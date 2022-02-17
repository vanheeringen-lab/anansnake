from snakemake.io import expand


rule motif2factors:
    """
    Create a gimme pfm/motif2factor file.
    For non-human/mouse, get ortholog TFs.
    """
    input:
        genome=GENOME,
    output:
        expand("{result_dir}/gimme/{assembly}.gimme.vertebrate.v5.0.pfm", assembly=ASSEMBLY, **config),
    log:
        expand("{result_dir}/gimme/log_{assembly}_m2f.txt", assembly=ASSEMBLY, **config),
    params:
        genomes_dir=config.get("genomes_dir"),
        tmpdir=config.get("tmp_dir"),
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
