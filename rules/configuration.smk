import os
import sys
from time import sleep

import pandas as pd
import genomepy
from snakemake.logging import logger


# all files found?
for item in ["samples", "gene_counts", "atac_counts"]:
    if item not in [None, "", "None"]:
        path = config[item]
        if not os.path.exists(path):
            logger.error(f"could not find '{item}' in {path}")
            sys.exit(1)
genomepy.Genome(config["genome"])

# read samples (contains column names & conditions)
samples = pd.read_csv(config["samples"], sep='\t', dtype='str', comment='#')

# error checking
for condition in samples.condition.unique():
    condition_samples = samples[samples["condition"] == condition]["sample"]

    # at least 2 samples per condition (for DESeq2)
    l = len(set(condition_samples))
    if l < 2:
        logger.error(f"condition {condition} had fewer than the minimum 2 samples")
        sys.exit(1)

    # no duplicate samples per condition
    if any(condition_samples.duplicated()):
        logger.error(f"duplicate samples in condition {condition}")
        sys.exit(1)

# make paths absolute
for item in ["result_dir", "gene_counts", "atac_counts"]:
    if item not in [None, "None"]:
        config[item] = genomepy.utils.cleanpath(config[item])
genomepy.utils.mkdir_p(config["result_dir"])

# print config
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")

sleep(1.5)
