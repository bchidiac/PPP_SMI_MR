#!/usr/bin/env python
# ---- Step 1a | UKB-PPP cis extraction ---------------------------------------
# Downloads each UKB-PPP protein GWAS from Synapse and reduces it to its cis
# variants (Supplementary Methods SM1, SM2):
#   - window: gene start - 1 Mb to gene end + 1 Mb on GRCh38, from the Olink
#     gene map, taken as the union over the protein's genes
#   - filters: INFO > 0.8 and 0.05 < A1FREQ < 0.95; p taken as 10^-LOG10P
#   - chromosome X is relabelled "X" to match the Olink map
#
# Reads   a Synapse folder (discovery or combined), Olink protein-to-gene map
# Writes  one cis file per protein; each archive is deleted after reduction
# Next    1a_join_exposure_rsid.py
# Needs   a Synapse account (~/.synapseConfig); on a cluster run as an array
#         job with hpc/1a_extract_exposure.sh
import argparse, os, re, sys, tarfile, time
import pandas as pd
import synapseclient

FOLDERS = {"discovery": "syn51365303", "combined": "syn51365308"}
COLS = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "LOG10P"]

p = argparse.ArgumentParser()
p.add_argument("--release", choices=FOLDERS, required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--task", type=int, default=0)
p.add_argument("--ntasks", type=int, default=1)
p.add_argument("--names-file", help="restrict the run to the archive names listed in this file")
p.add_argument("--gene-map", default="../MR_data/olink_protein_map_3k_v1.tsv")
a = p.parse_args()

gene_map = pd.read_csv(a.gene_map, sep="\t", dtype=str)
# the Olink map drops a digit from this gene_end (16419771), which would empty the cis window
gene_map.loc[gene_map.Assay == "ANP32C", "gene_end"] = "164197711"
gene_map["gene_start"] = gene_map.gene_start.astype(int)
gene_map["gene_end"] = gene_map.gene_end.astype(int)
case2_ids = set(gene_map.loc[gene_map.Assay.str.contains("_"), "OlinkID"])

clean_dir = os.path.join(a.outdir, "cis_clean")
tmp_dir = os.path.join(a.outdir, "tmp", f"task{a.task}")
os.makedirs(clean_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)

syn = synapseclient.Synapse(silent=True)
syn.login(silent=True)
kids = sorted(syn.getChildren(FOLDERS[a.release]), key=lambda k: k["name"])
if a.names_file:
    keep = set(open(a.names_file).read().split())
    kids = [k for k in kids if k["name"] in keep]
kids = kids[a.task::a.ntasks]
print(f"task {a.task}/{a.ntasks}: {len(kids)} archives", flush=True)

for i, k in enumerate(kids):
    stem = k["name"][:-4]
    out = os.path.join(clean_dir, stem + ".csv.gz")
    if os.path.exists(out):
        continue
    t0 = time.time()
    oid = re.search(r"OID\d+", stem).group(0)
    rows = gene_map[gene_map.OlinkID == oid]
    if rows.empty:
        print(f"[{i}] no map rows for {oid}; skipped", flush=True)
        continue
    if oid in case2_ids:
        gene, protein = rows.Assay.iloc[0], rows.UniProt.iloc[0]
    else:
        gene, protein = stem.split("_")[0], stem.split("_")[1]

    for attempt in range(5):
        try:
            ent = syn.get(k["id"], downloadLocation=tmp_dir, ifcollision="overwrite.local")
            break
        except Exception as err:
            print(f"[{i}] download attempt {attempt + 1} failed: {err}", flush=True)
            time.sleep(60 * (attempt + 1))
    else:
        sys.exit(f"[{i}] download failed for {oid}")

    parts = []
    with tarfile.open(ent.path) as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        for chrom, r in rows.groupby("chr", sort=False):
            prefix = f"{a.release}_chr{chrom}_"
            hit = [m for m in members if os.path.basename(m.name).startswith(prefix)]
            if len(hit) != 1:
                sys.exit(f"[{i}] {oid}: expected 1 member for chr{chrom}, found {len(hit)}")
            df = pd.read_csv(tf.extractfile(hit[0]), sep=" ", compression="gzip",
                             usecols=COLS, dtype={"CHROM": str, "ID": str})
            inwin = pd.Series(False, index=df.index)
            for s, e in zip(r.gene_start, r.gene_end):
                inwin |= (df.GENPOS >= s - 1_000_000) & (df.GENPOS <= e + 1_000_000)
            df = df[inwin].copy()
            df["CHROM"] = str(chrom)
            parts.append(df)
    os.remove(ent.path)

    cis = pd.concat(parts, ignore_index=True)
    n_cis = len(cis)
    cis = cis[(cis.INFO > 0.8) & (cis.A1FREQ > 0.05) & (cis.A1FREQ < 0.95)].copy()
    cis["P"] = 10 ** (-cis.LOG10P)
    cis = cis.drop(columns=["LOG10P", "GENPOS"]).rename(columns={"CHROM": "seqnames"})
    cis["dataset"] = stem + ".csv.gz"
    cis["gene"] = gene
    cis["protein"] = protein
    cis.to_csv(out + ".part", index=False, compression="gzip")
    os.replace(out + ".part", out)
    print(f"[{i}] {oid}: cis rows {n_cis}, kept {len(cis)}, {time.time() - t0:.0f}s", flush=True)

print(f"task {a.task} done", flush=True)
