# ---- Step 6 | Druggability --------------------------------------------------
# Annotates the MR-significant proteins from downloaded public tables, matched
# offline (Methods, Druggability):
#   - Open Targets tractability, grouped as in Dardani et al. (group 1, an
#     approved or clinical-stage drug; group 2, likely or potentially
#     druggable; group 3, little or unknown druggability)
#   - drugs with a ChEMBL mechanism of action on the protein, by highest
#     clinical stage
#   - Human Protein Atlas location: the secretome classification where a
#     protein is secreted, the predicted protein class otherwise
#   - DGIdb checked for approved drugs Open Targets does not list, written to
#     a separate check file
#
# Reads   mr_results_all_* from step 2, the Open Targets, HPA and DGIdb
#         downloads, Olink gene map
# Writes  one row per significant protein, and the DGIdb check
# Feeds   7 (Supplementary Table 5)

suppressMessages({library(data.table); library(arrow)})
options(arrow.skip_nul = TRUE)

# ---- Settings ----
env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
data_dir    <- env("MR_DATA_DIR", "../MR_data")
drug_dir    <- env("DRUG_DIR", file.path(data_dir, "druggability"))
ot_dir      <- env("OT_DIR", file.path(drug_dir, "OpenTargets_26.06"))
hpa_file    <- env("HPA_FILE", file.path(drug_dir, "HPA", "proteinatlas.tsv.zip"))
dgidb_file  <- env("DGIDB_FILE", file.path(drug_dir, "DGIdb", "interactions.tsv"))
olink_map   <- env("OLINK_MAP", file.path(data_dir, "olink_protein_map_3k_v1.tsv"))
results_dir <- env("MR_OUT_DIR", "results")
run_date    <- env("MR_RUN_DATE", as.character(Sys.Date()))

# ---- Proteins ----
files <- list.files(results_dir, paste0("^mr_results_all_.*_", run_date, "\\.csv$"), full.names = TRUE)
if (!length(files)) stop("No mr_results_all_* tables for ", run_date, " in ", results_dir)
sig <- unique(rbindlist(lapply(files, fread, select = c("id.exposure", "significant")))[significant == TRUE, .(protein = id.exposure)])

map <- unique(fread(olink_map)[, .(protein = sub(":.*", "", UKBPPP_ProteinID), hgnc = HGNC.symbol,
                                   ensembl_id, uniprot = UniProt)], by = "protein")
prot <- merge(sig, map, by = "protein", all.x = TRUE)

# ---- Open Targets tractability ----
# group 1: approved/clinical drug (any modality); 2: small-molecule buckets 4-8
# or antibody buckets 4-5; 3: otherwise
clinical_buckets <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
sm_group2 <- c("Structure with Ligand", "High-Quality Ligand", "High-Quality Pocket",
               "Med-Quality Pocket", "Druggable Family")                      # SM buckets 4-8
ab_group2 <- c("UniProt loc high conf", "GO CC high conf")                    # AB buckets 4-5
tgt <- as.data.table(open_dataset(file.path(ot_dir, "target")) |>
                       dplyr::select(id, tractability) |> dplyr::collect())
tgt <- tgt[id %in% prot$ensembl_id]
tract <- rbindlist(lapply(seq_len(nrow(tgt)), function(k) {
  x <- tgt$tractability[[k]]
  if (!length(x)) return(NULL)
  data.table(ensembl_id = tgt$id[k], modality = x$modality, bucket = x$id, value = x$value)
}))[value == TRUE]
best <- function(m) tract[modality == m, .(label = bucket[1]), by = ensembl_id]   # buckets listed in rank order
tract_dt <- tract[, .(
  tractability_group = if (any(bucket %in% clinical_buckets)) 1L
                       else if (any((modality == "SM" & bucket %in% sm_group2) |
                                    (modality == "AB" & bucket %in% ab_group2))) 2L else 3L),
  by = ensembl_id]
for (m in c("SM", "AB", "PR", "OC"))
  tract_dt <- merge(tract_dt, setnames(best(m), "label", paste0("best_", tolower(m))),
                    by = "ensembl_id", all.x = TRUE)

# ---- Drugs acting on each protein, by highest clinical stage ----
stage_rank <- c(APPROVAL = 10, PREAPPROVAL = 9, PHASE_3 = 8, PHASE_2_3 = 7, PHASE_2 = 6,
                PHASE_1_2 = 5, PHASE_1 = 4, EARLY_PHASE_1 = 3, IND = 2, PRECLINICAL = 1, UNKNOWN = 0)
ct <- as.data.table(open_dataset(file.path(ot_dir, "clinical_target")) |>
                      dplyr::select(drugId, targetId, maxClinicalStage) |> dplyr::collect())
ct <- ct[targetId %in% prot$ensembl_id]
mol <- as.data.table(open_dataset(file.path(ot_dir, "drug_molecule")) |>
                       dplyr::select(id, name) |> dplyr::collect())
ct <- merge(ct, mol, by.x = "drugId", by.y = "id", all.x = TRUE)
ct[, rank := stage_rank[maxClinicalStage]]
drugs_dt <- ct[, .(
  max_clinical_stage = maxClinicalStage[which.max(rank)],
  n_drugs_approved   = uniqueN(drugId[maxClinicalStage == "APPROVAL"]),
  n_drugs_clinical   = uniqueN(drugId[rank >= 3 & maxClinicalStage != "APPROVAL"]),
  approved_drugs     = paste(sort(unique(name[maxClinicalStage == "APPROVAL"])), collapse = "; ")),
  by = .(ensembl_id = targetId)]

# ---- Human Protein Atlas location ----
# HPA's secretome field lumps every non-secreted protein into one category, hence the protein class
# for those.
hpa <- fread(cmd = paste("unzip -p", shQuote(hpa_file)),
             select = c("Ensembl", "Secretome location", "Protein class"))
setnames(hpa, c("ensembl_id", "secretome", "protein_class"))
pred <- function(x) vapply(strsplit(hpa$protein_class, ","), function(p)
  any(trimws(p) == paste0("Predicted ", x, " proteins")), logical(1))
hpa[, hpa_location := fcase(
  grepl("^Secreted", secretome), secretome,
  pred("intracellular") & pred("membrane"), "Intracellular and membrane",
  pred("intracellular"), "Intracellular",
  pred("membrane"), "Membrane",
  pred("secreted"), "Secreted",
  default = NA_character_)]
hpa <- hpa[, .(ensembl_id, hpa_location)]

# ---- DGIdb check ----
# Included only sources that record drug-target pharmacology; the rest record
# cancer variant-therapy associations. Reported results remained the same
# without this restriction.
curated <- c("ChEMBL", "GuideToPharmacology", "TTD", "PharmGKB", "FDA", "TEND", "TdgClinicalTrial")
dg <- fread(dgidb_file, quote = "")[interaction_source_db_name %in% curated & approved == TRUE &
                                     gene_name %in% c(prot$protein, prot$hgnc)]
dg <- merge(dg, unique(rbind(prot[, .(gene_name = protein, ensembl_id)], prot[, .(gene_name = hgnc, ensembl_id)])),
            by = "gene_name")
dg <- dg[, .(sources = paste(sort(unique(interaction_source_db_name)), collapse = "; "),
             pharmgkb_only = all(interaction_source_db_name == "PharmGKB")),
         by = .(ensembl_id, drug_uc = toupper(drug_name))]
dg <- dg[!paste(ensembl_id, drug_uc) %in% ct[, paste(targetId, toupper(name))]]   # not listed by Open Targets
mol_all <- as.data.table(open_dataset(file.path(ot_dir, "drug_molecule")) |>
                           dplyr::select(id, name, synonyms, tradeNames) |> dplyr::collect())
names_lut <- unique(rbind(
  mol_all[, .(id, drug_uc = toupper(name))],
  mol_all[, .(drug_uc = toupper(unlist(lapply(synonyms, function(x) x$label)))), by = id],
  mol_all[, .(drug_uc = toupper(unlist(lapply(tradeNames, function(x) x$label)))), by = id])[!is.na(drug_uc)],
  by = "drug_uc")
moa <- as.data.table(open_dataset(file.path(ot_dir, "drug_mechanism_of_action")) |>
                       dplyr::select(chemblIds, targets) |> dplyr::collect())
moa <- moa[, .(id = unlist(chemblIds), target = rep(list(unlist(targets)), length(unlist(chemblIds)))), by = seq_len(nrow(moa))]
moa <- moa[, .(target = unlist(target)), by = id]
dg <- merge(dg, names_lut, by = "drug_uc", all.x = TRUE)
dg[, category := fifelse(pharmgkb_only, "pharmacogenomic association only (PharmGKB)",
                 fifelse(is.na(id), "drug not in ChEMBL",
                 fifelse(paste(id, ensembl_id) %in% moa[, paste(id, target)], "mechanism of action on this protein",
                 fifelse(id %in% moa$id, "mechanism of action on a different protein",
                         "no mechanism of action recorded"))))]
dgidb_check <- merge(dg[, .(ensembl_id, drug = drug_uc, sources, category)], prot[, .(ensembl_id, protein)], by = "ensembl_id")
fwrite(dgidb_check[order(protein, category, drug)],
       file.path(results_dir, paste0("druggability_dgidb_check_", run_date, ".csv")))

# ---- Output ----
out <- Reduce(function(a, b) merge(a, b, by = "ensembl_id", all.x = TRUE),
              list(prot, tract_dt, drugs_dt, hpa))
out[is.na(tractability_group) & !is.na(ensembl_id), tractability_group := 3L]
out[is.na(n_drugs_approved), `:=`(n_drugs_approved = 0L, n_drugs_clinical = 0L)]
setcolorder(out, c("protein", "hgnc", "uniprot", "ensembl_id"))
fwrite(out[order(protein)], file.path(results_dir, paste0("druggability_", run_date, ".csv")))

cat("Proteins annotated:", nrow(out), "; not matched to Ensembl:", sum(is.na(out$ensembl_id)), "\n")
cat("Tractability group 1 / 2 / 3:", paste(tabulate(out$tractability_group, 3), collapse = " / "), "\n")
cat("With an approved drug:", sum(out$n_drugs_approved > 0), "; with a clinical-stage drug only:",
    sum(out$n_drugs_approved == 0 & out$n_drugs_clinical > 0), "\n")
cat("HPA location:\n"); print(out[, .N, by = hpa_location][order(-N)])
cat("DGIdb check - approved drugs not listed by Open Targets, by category:\n"); print(dgidb_check[, .N, by = category])
cat("Proteins with a missed drug whose mechanism of action is the protein:",
    uniqueN(dgidb_check[category == "mechanism of action on this protein", protein]), "\n")
