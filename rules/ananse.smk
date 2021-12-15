from snakemake.io import expand


rule binding:
    """
    ANANSE binding
    """
    input:
        atac=config["atac_counts"],
        pfm=rules.motif2factor.output,
        pfmscorefile=rules.pfmscorefile.output,
    output:
        expand("{result_dir}/binding/binding.h5", **config),
    log:
        expand("{result_dir}/binding/log.txt", **config),
    benchmark:
        expand("{result_dir}/benchmarks/binding.txt", **config)[0]
    params:
        genome=config["genome"],
        jaccard=config["jaccard"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})
        
        ananse binding \
        -A {input.atac} \
        -g {params.genome} \
        -p {input.pfm} \
        --pfmscorefile {input.pfmscorefile} \
        --jaccard-cutoff {params.jaccard} \
        -n {threads} \
        -o (dirname {output}) \
        > {log} 2>&1
        """

# PARAMS
# hist = config["hist_counts"],
# regions = config["regions"],
# reference = config["reference"],
# samples=config["binding_samples"],

# SHELL
# -H {params.hist} \
# -r {params.regions} \
# -R {params.reference} \
# -c {params.samples} \

rule network:
    """
    ANANSE network
    """
    input:
        binding=rules.binding.output,
        genes=config["gene_counts"],
    output:
        expand("{result_dir}/network/{{condition}}.tsv",**config),
    log:
        expand("{result_dir}/network/log_{{condition}}.txt",**config),
    benchmark:
        expand("{result_dir}/benchmarks/network_{{condition}}.txt",**config)[0]
    params:
        samples=lambda wildcards: list(samples[samples["condition"] == wildcards.condition]["sample"]),
        genome=config["genome"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})

        ananse network \
        {input.binding} \
        -e {input.genes} \
        -c {params.samples} \
        -g {params.genome} \
        -o {output} \
        --full-output \
        -n {threads} \
        > {log} 2>&1
        """

def get_contrasts(wildcards):
    networks = dict()
    column, target, source = wildcards.contrast.split("_")
    networks["target"] = f"{config['result_dir']}/network/{target}.tsv"
    networks["source"] = f"{config['result_dir']}/network/{source}.tsv"
    return networks

rule influence:
    """
    ANANSE influence
    """
    input:
        unpack(get_contrasts),
        degenes=rules.deseq2.output,
    output:
        inf = expand("{result_dir}/influence/{{contrast}}.tsv",**config),
        diff_inf = expand("{result_dir}/influence/{{contrast}}_diffnetwork.tsv",** config),
    log:
        expand("{result_dir}/influence/log_{{contrast}}.txt",**config),
    benchmark:
        expand("{result_dir}/benchmarks/influence_{{contrast}}.txt",**config)[0]
    params:
        genome=config["genome"],
        edges=config["edges"],
        padj=config["padj"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})

        ananse influence \
        -s {input.source} \
        -t {input.target} \
        -d {input.degenes} \
        -a {params.genome} \
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
    benchmark:
        expand("{result_dir}/benchmarks/plot_{{contrast}}.txt",**config)[0]
    params:
        type=config["plot_type"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        # for the log
        mkdir -p (dirname {output})

        ananse plot \
        {input.inf} \
        -d {input.diff_inf} \
        -t {params.type} \
        --full-output \
        -o (dirname {output}) \
        > {log} 2>&1
        """
