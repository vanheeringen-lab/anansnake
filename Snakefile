from snakemake.io import expand


include: "rules/configuration.smk"
include: "rules/gimme.smk"
include: "rules/deseq2.smk"
include: "rules/ananse.smk"


rule all:
    input:
        expand("{result_dir}/plot/{contrasts}/influence.{plot_type}", **config)
