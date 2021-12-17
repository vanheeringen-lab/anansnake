import os
import sys
from time import sleep

import pandas as pd
import genomepy
from snakemake.logging import logger
from seq2science.util import parse_samples, parse_contrast

# NOTE: global variables are written in all caps

# check + fix file paths
for item in ["result_dir", "rna_samples", "rna_counts", "rna_tpms", "atac_samples", "atac_counts"]:
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

# check if the genome is found + add to globals
g = genomepy.Genome(config["genome"], config.get("genomes_dir"), build_index=False)
GENOME = g.genome_dir
ASSEMBLY = g.name

# read samples (contains column names & conditions)
rna_samples = pd.read_csv(config["rna_samples"], sep='\t', dtype='str', comment='#')
atac_samples = pd.read_csv(config["atac_samples"], sep='\t', dtype='str', comment='#')
# mimic s2s samples files
# TODO: how to deal with unmerged replicates (technical/biological/merged in only 1 workflow)
rna_samples = parse_samples(rna_samples, config)
atac_samples = parse_samples(atac_samples, config)

CONDITIONS = dict()
CONTRASTS = dict()
sampledict = {
    "RNA-seq": rna_samples,
    "ATAC-seq": atac_samples,
}
for wf, samples in sampledict.items():
    # check for minimal columns
    assert "sample" in samples.columns
    if wf == "RNA-seq":
        assert "assembly" in samples.columns, "deseq2science requires an assembly column"

    # which column contains names present in the counts/TPM/peaks tables
    samplecol = "sample"
    if "technical_replicates" in samples:
        samplecol = "technical_replicates"
    elif "_trep" in samples:
        samplecol = "_trep"

    for contrast in config["contrasts"]:
        batch, col, target, source = parse_contrast(contrast, samples, check=True)
        # conditions per contrast (order is important)
        if wf == "RNA-seq":
            CONTRASTS[contrast] = {"target": target, "source": source}

        # samples & column per condition
        for condition in [target, source]:
            if condition not in CONDITIONS:
                CONDITIONS[condition] = {"column": [], "RNA-seq samples": [], "ATAC-seq samples": []}

            # the columns should overlap between samples files
            if col not in samples.columns:
                logger.error(f"Column '{col}' not found in the {wf} samples file!")
                sys.exit(1)

            if col not in CONDITIONS[condition]["column"]:
                CONDITIONS[condition]["column"].append(col)
            for sample in list(samples[samples[col] == condition][samplecol]):
                if sample not in CONDITIONS[condition][f"{wf} samples"]:
                    CONDITIONS[condition][f"{wf} samples"].append(sample)

# check if each condition occurs in only one used column (cant distinguish conditions between columns)
for condition in CONDITIONS:
    l = len(CONDITIONS[condition]["column"])
    if l == 0:
        # this really shouldn't happen
        logger.error(
            f"Contrast condition '{condition}' not found in any column."
        )
        sys.exit(1)
    if "_trep" in CONDITIONS[condition]["column"]:
        l -= 1
    if "_brep" in CONDITIONS[condition]["column"]:
        l -= 1
    if l > 1:
        logger.error(
            "Contrast conditions may not be used in multiple columns! ",
            f"'{condition}' occurs in more than one column ({CONDITIONS[condition]['column']})."
        )
        sys.exit(1)

    # minimum number of samples per condition
    if len(CONDITIONS[condition]["RNA-seq samples"]) < 2:
        logger.error(
            f"In contrast condition '{condition}', "
            "one of the conditions has fewer than the minimum 2 RNA-seq samples"
        )
        sys.exit(1)
    if len(CONDITIONS[condition]["ATAC-seq samples"]) < 1:
        logger.error(
            f"In contrast condition '{condition}', "
            "one of the conditions has no ATAC-seq samples"
        )
        sys.exit(1)

# print config
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")
sleep(1.5)

# set default resource limits
for res in ["deseq2", "network"]:
    if res not in workflow.global_resources:
        workflow.global_resources[res] = 1
