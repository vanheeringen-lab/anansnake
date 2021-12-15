import os
import sys
from time import sleep

import pandas as pd
import genomepy
from snakemake.logging import logger
from seq2science.util import parse_samples

# check + fix file paths
for item in ["result_dir", "rna_samples", "rna_counts", "atac_samples", "atac_counts"]:
    path = genomepy.utils.cleanpath(config[item])

    # create output dir
    if item == "result_dir":
        genomepy.utils.mkdir_p(path)

    # all files found?
    elif not os.path.exists(path):
        logger.error(f"could not find '{item}' in {path}")
        sys.exit(1)

    # save absolute path
    config[item] = path

# check if the genome is found
g = genomepy.Genome(config["genome"], build_index=False)
config["assembly"] = g.name

# read samples (contains column names & conditions)
rna_samples = pd.read_csv(config["rna_samples"], sep='\t', dtype='str', comment='#')
atac_samples = pd.read_csv(config["atac_samples"], sep='\t', dtype='str', comment='#')

# check columns
assert "sample" in rna_samples
assert "sample" in atac_samples
for contrast in config["contrasts"]:
    col, target, source = contrast.split("_")
    if col not in rna_samples:
        logger.error(f"For contrast '{contrast}', column '{col}' wasn't found in the RNA-seq samples '{rna_samples}'")
        sys.exit(1)
    if col not in atac_samples:
        logger.error(f"For contrast '{contrast}', column '{col}' wasn't found in the ATAC-seq samples '{atac_samples}'")
        sys.exit(1)

    # check rows
    tgt_samples = rna_samples[rna_samples[col] == target]["sample"]
    src_samples = rna_samples[rna_samples[col] == source]["sample"]
    for samples in [tgt_samples, src_samples]:
        # at least 2 samples per condition (for DESeq2)
        l = len(set(samples))
        if l < 2:
            logger.error(f"In contrast {contrast}, one of the conditions has fewer than the minimum 2 RNA-seq samples")
            sys.exit(1)

        # no duplicate samples
        if any(samples.duplicated()):
            logger.error(f"duplicate samples in contrast {contrast} in RNA-seq samples")
            sys.exit(1)

    tgt_samples = atac_samples[atac_samples[col] == target]["sample"]
    src_samples = atac_samples[atac_samples[col] == source]["sample"]
    for samples in [tgt_samples, src_samples]:
        # at least 1 sample per condition
        l = len(set(samples))
        if l == 0:
            logger.error(f"In contrast {contrast}, one of the conditions has no ATAC-seq samples")
            sys.exit(1)

        # no duplicate samples
        if any(samples.duplicated()):
            logger.error(f"duplicate samples in contrast {contrast} in ATAC-seq samples")
            sys.exit(1)

# mimic s2s samples files
if "assembly" not in rna_samples:
    rna_samples["assembly"] = g.name
rna_samples = parse_samples(rna_samples, config)
if "assembly" not in atac_samples:
    atac_samples["assembly"] = g.name
atac_samples = parse_samples(atac_samples, config)

# print config
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")
sleep(1.5)
