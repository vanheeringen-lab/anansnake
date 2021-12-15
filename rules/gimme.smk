from snakemake.io import expand


rule motif2factor:
    """
    Create a gimme pfm/motif2factor with ortholog TFs for this genome
    """
    output:
        expand("{result_dir}/gimme/{assembly}.gimme.vertebrate.v5.0.pfm",** config),
    log:
        expand("{result_dir}/gimme/log_{assembly}_m2f.txt",**config),
    params:
        genome=config["genome"],
    threads: 12
    conda: "../envs/gimme.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        gimme motif2factors \
        --new-reference {params.genome} \
        --outdir $outdir \
        --tmpdir {resources.tmpdir} \
        --threads {threads} \
        > {log} 2>&1
        """


rule pfmscorefile:
    """
    Scan motif activity
    - in all putative enhancer regions
    - for all (ortholog) TFs
    """
    input:
        regions=config["atac_counts"],
        pfm=rules.motif2factor.output
    output:
        expand("{result_dir}/gimme/pfmscorefile.tsv",**config),
    log:
        expand("{result_dir}/gimme/log_{assembly}_pfmscorefile.txt",**config),
    params:
        genome=config["genome"],
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
        -g {params.genome} \
        -Tz --gc \
        -n {threads} \
        > {output} \
        2> {log}
        """
