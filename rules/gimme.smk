from snakemake.io import expand


rule motif2factor:
    """
    gimme motif2factors
    """
    output:
        expand("{result_dir}/gimme/{genome}.gimme.vertebrate.v5.0.pfm",** config),
    log:
        expand("{result_dir}/gimme/log_{genome}_m2f.txt",**config),
    params:
        genome=config["genome"],
    threads: 12
    conda: "../envs/gimme.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})

        gimme motif2factors \
        --new-reference {params.genome} \
        --outdir (dirname {output}) \
        --tmpdir {resources.tmpdir} \
        --threads {threads} \
        > {log} 2>&1
        """


rule pfmscorefile:
    """
    gimme scam
    """
    input:
        regions=config["atac_counts"],
        pfm=rules.motif2factor.output
    output:
        expand("{result_dir}/gimme/pfmscorefile.tsv",**config),
    log:
        expand("{result_dir}/gimme/log_{genome}_pfmscorefile.txt",**config),
    params:
        genome=config["genome"],
    threads: 12
    conda: "../envs/gimme.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})

        gimme scan \
        {input.regions} \
        -p {input.pfm} \
        -g {params.genome} \
        -Tz --gc \
        -n {threads} \
        > {output} \
        2> {log}
        """
