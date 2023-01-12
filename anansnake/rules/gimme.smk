from snakemake.io import expand


rule motif2factors:
    """
    Create a gimme pfm/motif2factor file.
    If get_orthologs is True, find ortholog TFs, else copy the database files as-is.
    """
    input:
        genome=GENOME,
    output:
        expand("{result_dir}/gimme/{assembly}.{database}.pfm", assembly=ASSEMBLY, **config),
    log:
        expand("{result_dir}/gimme/log_{assembly}_m2f.txt", assembly=ASSEMBLY, **config),
    params:
        genomes_dir=config.get("genomes_dir"),
        database=config.get("database", "gimme.vertebrate.v5.0"),
        get_orthologs=config.get("get_orthologs", True),
        tmpdir=config.get("tmp_dir"),
        keeptmp=config.get("keep_tmp_data", False),
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
        expand("{result_dir}/gimme/pfmscorefile.tsv", **config),
    log:
        expand("{result_dir}/gimme/log_{assembly}_pfmscorefile.txt", assembly=ASSEMBLY, **config),
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
        directory(expand("{result_dir}/gimme/{assembly}-maelstrom", assembly=ASSEMBLY, **config)),
    log:
        expand("{result_dir}/gimme/log_{assembly}_maelstrom.txt", assembly=ASSEMBLY, **config),
    params:
        atac_samples=lambda wildcards : sorted({sample for vals in CONDITIONS.values() for sample in vals}),
    threads: 24
    resources:
        mem_mb=40_000,
    conda: "../envs/gimme.yaml"
    script:
        "../scripts/maelstrom.py"
