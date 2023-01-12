from snakemake.io import expand


rule binding:
    """
    Measure enhancer activity for one specific condition
    """
    input:
        atac=config["atac_counts"],
        pfm=rules.motif2factors.output,
        pfmscorefile=rules.pfmscorefile.output,
        genome=GENOME,
    output:
        expand("{result_dir}/binding/{{condition}}.h5", ** config),
    log:
        expand("{result_dir}/binding/log_{{condition}}.txt", **config),
    benchmark:
        expand("{result_dir}/benchmarks/binding_{{condition}}.txt", **config)[0]
    params:
        atac_samples=lambda wildcards: CONDITIONS[wildcards.condition]["ATAC-seq samples"],
        enhancer_data=lambda wildcards: "-P" if config.get("enhancer_data") == "p300" else "-A",
        jaccard=config["jaccard"],
    threads: 1  # multithreading not required when using a pfmscorefile
    resources:
        mem_mb=40_000,  # 30-50 GB
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})/{wildcards.condition}
        trap "rm -rf $outdir;" EXIT

        # for the log
        mkdir -p $outdir

        # additional log info
        printf "using columns: {params.atac_samples}\n\n" > {log}

        ananse binding \
        {params.enhancer_data} {input.atac} \
        -c {params.atac_samples} \
        -g {input.genome} \
        -p {input.pfm} \
        --pfmscorefile {input.pfmscorefile} \
        --jaccard-cutoff {params.jaccard} \
        -n {threads} \
        -o $outdir \
        >> {log} 2>&1
        
        if [ -f $outdir/binding.h5 ]; then
            mv $outdir/binding.h5 {output}
        fi
        """

# PARAMS
# hist = config["hist_counts"],
# regions = config["regions"],
# reference = config["reference"],

# SHELL
# -H {params.hist} \
# -r {params.regions} \
# -R {params.reference} \


rule network:
    """
    generate a GRN for one specific condition
    """
    input:
        binding=rules.binding.output,
        genes=config["rna_tpms"],
        genome=GENOME,
    output:
        expand("{result_dir}/network/{{condition}}.tsv",**config),
    log:
        expand("{result_dir}/network/log_{{condition}}.txt",**config),
    benchmark:
        expand("{result_dir}/benchmarks/network_{{condition}}.txt",**config)[0]
    params:
        rna_samples=lambda wildcards: CONDITIONS[wildcards.condition]["RNA-seq samples"],
    threads: 1  # multithreading explodes memory
    resources:
        network=1,
        mem_mb=24_000,  # 8-24 GB
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        # additional log info
        printf "using columns: {params.rna_samples}\n\n" > {log}

        ananse network \
        {input.binding} \
        -e {input.genes} \
        -c {params.rna_samples} \
        -g {input.genome} \
        -o {output} \
        --full-output \
        -n {threads} \
        >> {log} 2>&1
        """


def get_source(wildcards):
    return f"{config['result_dir']}/network/{CONTRASTS[wildcards.contrast]['source']}.tsv"

def get_target(wildcards):
    return f"{config['result_dir']}/network/{CONTRASTS[wildcards.contrast]['target']}.tsv"

rule influence:
    """
    Find the most influential TFs between to conditions
    - as defined in the config.yaml, under contrasts.
    """
    input:
        source=get_source,
        target=get_target,
        degenes=rules.deseq2.output,
        genome=GENOME,
    output:
        inf = expand("{result_dir}/influence/{{contrast}}.tsv",**config),
        diff_inf = expand("{result_dir}/influence/{{contrast}}_diffnetwork.tsv",** config),
    log:
        expand("{result_dir}/influence/log_{{contrast}}.txt",**config),
    benchmark:
        expand("{result_dir}/benchmarks/influence_{{contrast}}.txt",**config)[0]
    params:
        edges=config["edges"],
        padj=config["padj"],
    threads: 1  # multithreading not required with the new networkx implementation
    resources:
        mem_mb=24_000,  # 7-150 GB
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        ananse influence \
        -s {input.source} \
        -t {input.target} \
        -d {input.degenes} \
        -a {input.genome} \
        -i {params.edges} \
        -j {params.padj} \
        --full-output \
        -o {output.inf} \
        -n {threads} \
        > {log} 2>&1
        """


rule plot:
    """
    ANANSE plot
    """
    input:
        inf=rules.influence.output.inf,
        diff_inf=rules.influence.output.diff_inf,
    output:
        expand("{result_dir}/plot/{{contrast}}.{plot_type}",**config),
    log:
        expand("{result_dir}/plot/log_{{contrast}}.txt",**config),
    params:
        type=config["plot_type"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})/{wildcards.contrast}
        trap "rm -rf $outdir;" EXIT

        # for the log
        mkdir -p $outdir

        ananse plot \
        {input.inf} \
        -d {input.diff_inf} \
        -t {params.type} \
        --full-output \
        -o $outdir \
        > {log} 2>&1
        
        if [ -f $outdir/influence.{params.type} ]; then
            mv $outdir/influence.{params.type} {output}
        fi
        if [ -f $outdir/topTF_GRN.{params.type} ]; then
            mv $outdir/topTF_GRN.{params.type} $(dirname {output})/{wildcards.contrast}_topTF_GRN.{params.type}
        fi
        """
