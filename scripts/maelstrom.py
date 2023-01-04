from contextlib import redirect_stdout, redirect_stderr
import os
import logging

import numpy as np
import pandas as pd
import qnorm

from gimmemotifs.maelstrom import run_maelstrom


log = snakemake.log
regions = snakemake.input.regions
genome = snakemake.input.genome
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
        # remove average column (if present)
        df = pd.read_table(regions, sep="\t", comment="#", index_col=0)
        if "average" in df.columns:
            df.drop(columns=["average"], inplace=True)

        # normalize count matrix
        df = np.log2(df + 1)
        df = qnorm.quantile_normalize(df)

        #reduce df lenght to only highly variable rows if df length > 50.000
        if len(df) > 50000:
            df['variance'] = df.var(axis = 1)
            df = df.sort_values('variance', ascending = False)[1:50000]
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
