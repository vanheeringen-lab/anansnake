from contextlib import redirect_stdout, redirect_stderr
import os
import logging
import re

import numpy as np
import pandas as pd
import qnorm

from gimmemotifs.maelstrom import run_maelstrom


log = snakemake.log
regions = snakemake.input.regions
genome = snakemake.input.genome
columns = snakemake.params.atac_samples
outdir = snakemake.output[0]
pfmfile = snakemake.input.pfm[0]
ncpus = snakemake.threads


# log gimme messages
logger = logging.getLogger("gimme")
logger.handlers.clear()
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
handler = logging.FileHandler(str(log))
handler.setLevel(logging.INFO)
handler.setFormatter(formatter)
logger.addHandler(handler)

# log other messages
with open(str(log), "a") as f:
    with redirect_stdout(f), redirect_stderr(f):
        df = pd.read_table(regions, sep="\t", comment="#", index_col=0)

        # AnanseSeurat/Scanpy: remove average column
        if "average" in df.columns:
            df.drop(columns=["average"], inplace=True)

        # filter for samples in the conditions
        logger.info(f"columns found in ATAC-seq read counts table: {df.columns}")
        cols = "|".join(columns)
        re_column = re.compile(rf"^{cols}$", re.IGNORECASE)
        df = df.filter(regex=re_column)
        logger.info(f"columns used for maelstrom: {df.columns}")
        if len(columns) != len(df.columns):
            logger.warning(
                f"{len(columns)} expected in ATAC-seq read counts table, "
                f"{len(df.columns)} found after filtering."
            )

        # remove zero rows (introduced by filtering columns)
        df = df[df.values.sum(axis=1) != 0]

        # normalize count matrix
        df = np.log2(df + 1)
        df = qnorm.quantile_normalize(df)

        # only use highly variable regions
        if len(df) > 50_000:
            df['variance'] = df.var(axis = 1)
            df = df.sort_values('variance', ascending = False).head(50_000)
            df = df.drop('variance', axis = 1)

        # save count matrix to temporary file
        os.makedirs(outdir, exist_ok=True)
        tmp = os.path.join(outdir, "maelstrom_input.tsv")
        df.to_csv(tmp, sep="\t")

        # gimme maelstrom
        run_maelstrom(
            infile=tmp,
            genome=genome,
            outdir=outdir,
            pfmfile=pfmfile,
            ncpus=ncpus,
            plot_all_motifs=True,
            random_state=np.random.RandomState(123),
        )
