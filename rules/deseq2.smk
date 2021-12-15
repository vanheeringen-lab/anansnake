from snakemake.io import expand


rule deseq2:
    """
    run DESeq2 using seq2science
    """
    input:
        rna_samples = config["rna_samples"],
        genes = config["rna_counts"],
    output:
        expand("{result_dir}/deseq2/{{contrast}}.diffexp.tsv",**config),
    log:
        expand("{result_dir}/deseq2/log_{{contrast}}.txt",**config),
    conda: "../envs/deseq2.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})
        
        deseq2science \
        -d {wildcards.contrast} \
        -s {input.rna_samples} \
        -c {input.genes} \
        -o (dirname {output}) \
        > {log} 2>&1
        """
