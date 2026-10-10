#!/usr/bin/env python
# ---- Step 1a | rsID join ----------------------------------------------------
# Combines the per-protein cis files and adds rsIDs and GRCh37/GRCh38 positions
# by an inner join on variant ID (Supplementary Methods SM1).
#
# Reads   the cis files from 1a_extract_exposure_cis.py, and the UKB-PPP rsID
#         map reduced to one record per variant ID, with columns ID, rsid,
#         POS19 and POS38 (read in chunks: it covers the whole genome)
# Writes  the exposure file read by step 2
import argparse, glob, os
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--indir", required=True)
p.add_argument("--out", required=True)
p.add_argument("--rsid-dict", default="../MR_data/ukbppp_rsid_maps/rsids_all_chromosomes_uniq2.tsv.gz")
a = p.parse_args()

files = sorted(glob.glob(os.path.join(a.indir, "*.csv.gz")))
allfiles = pd.concat((pd.read_csv(f, dtype={"seqnames": str, "ID": str}) for f in files), ignore_index=True)
print(f"files {len(files)}, rows before join {len(allfiles)}", flush=True)

ids = pd.Index(allfiles.ID.unique())
dict_parts = []
for chunk in pd.read_csv(a.rsid_dict, dtype={"ID": str}, usecols=["ID", "rsid", "POS19", "POS38"], chunksize=5_000_000):
    dict_parts.append(chunk[chunk.ID.isin(ids)])
rsid = pd.concat(dict_parts, ignore_index=True)

allfiles = allfiles.merge(rsid, on="ID", how="inner")
allfiles = allfiles.rename(columns={
    "seqnames": "chr.exposure", "POS19": "pos.exposure.19", "POS38": "pos.exposure.38", "rsid": "id.exposure",
    "ALLELE0": "other_allele.exposure", "ALLELE1": "effect_allele.exposure", "A1FREQ": "eaf.exposure",
    "N": "samplesize.exposure", "BETA": "beta.exposure", "SE": "se.exposure", "P": "pval.exposure"})
allfiles = allfiles[["chr.exposure", "pos.exposure.19", "pos.exposure.38", "id.exposure", "other_allele.exposure",
                     "effect_allele.exposure", "eaf.exposure", "INFO", "samplesize.exposure", "beta.exposure",
                     "se.exposure", "pval.exposure", "dataset", "gene", "protein"]]
print(f"rows after join {len(allfiles)}", flush=True)
allfiles.to_csv(a.out + ".part", sep="\t", index=False, compression="gzip")
os.replace(a.out + ".part", a.out)
