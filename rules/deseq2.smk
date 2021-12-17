from snakemake.io import expand


rule deseq2:
    """
    run DESeq2 using seq2science
    """
    input:
        rna_samples = config["rna_samples"],
        genes = config["rna_counts"],
    output:
        expand("{result_dir}/deseq2/{assembly}-{{contrast}}.diffexp.tsv", assembly=ASSEMBLY, **config),
    log:
        expand("{result_dir}/deseq2/log_{{contrast}}.txt", **config),
    resources:
        deseq2=1
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir
        
        deseq2science \
        -d {wildcards.contrast} \
        -s {input.rna_samples} \
        -c {input.genes} \
        -o $outdir \
        > {log} 2>&1
        """
