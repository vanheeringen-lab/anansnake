from snakemake.io import expand


rule motif2factors:
    """
    Create a gimme pfm/motif2factor file.
    If get_orthologs is True, find ortholog TFs, else copy the database files as-is.
    """
    input:
        genome=GENOME,
    output:
        expand("{gimme_dir}/{assembly}.{database}.pfm", assembly=ASSEMBLY, **config),
    log:
        expand("{log_dir}/gimme_m2f.txt", **config),
    params:
        genomes_dir=config.get("genomes_dir"),
        database=config.get("database", "gimme.vertebrate.v5.0"),
        get_orthologs=config.get("get_orthologs", True),
        keep_tmp=config.get("keep_ortholog_data", False),
    threads: 24
    conda: "../envs/gimme.yaml"
    script:
        "../scripts/motif2factors.py"


rule pfmscorefile:
    """
    Scan motif activity
    - in all putative enhancer regions
    - for all (ortholog) TFs
    """
    input:
        regions=config["atac_counts"],
        pfm=rules.motif2factors.output,
        genome=GENOME,
    output:
        expand("{gimme_dir}/pfmscorefile.tsv", **config),
    log:
        expand("{log_dir}/gimme_pfmscorefile.txt", **config),
    threads: 12
    conda: "../envs/gimme.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        gimme scan \
        {input.regions} \
        -p {input.pfm} \
        -g {input.genome} \
        -Tz --gc \
        -N {threads} \
        > {output} \
        2> {log}
        """


rule maelstrom:
    """
    Find differential motifs
    """
    input:
        regions=config["atac_counts"],
        pfm=rules.motif2factors.output,
        genome=GENOME,
    output:
        directory(expand("{gimme_dir}/{assembly}-maelstrom", assembly=ASSEMBLY, **config)),
    log:
        expand("{log_dir}/gimme_maelstrom.txt", **config),
    params:
        atac_samples=lambda wildcards : sorted({sample for v in CONDITIONS.values() for sample in v['ATAC-seq samples']}),
    threads: 24
    resources:
        mem_mb=40_000,
    conda: "../envs/gimme.yaml"
    script:
        "../scripts/maelstrom.py"
