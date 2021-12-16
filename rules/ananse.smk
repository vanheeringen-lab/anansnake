from snakemake.io import expand, unpack
from seq2science.util import parse_contrast


def get_samples(samples, contrasts, condition):
    ret = []
    cols = [parse_contrast(c, samples, check=False)[1] for c in contrasts]
    for col in cols:
        ret.extend(list(samples[samples[col] == condition]["sample"]))
    return list(set(ret))


rule binding:
    """
    Measure enhancer activity for one specific condition
    """
    input:
        atac=config["atac_counts"],
        pfm=rules.motif2factor.output,
        pfmscorefile=rules.pfmscorefile.output,
    output:
        expand("{result_dir}/binding/{{condition}}/binding.h5", **config),
    log:
        expand("{result_dir}/binding/log_{{condition}}.txt", **config),
    benchmark:
        expand("{result_dir}/benchmarks/binding_{{condition}}.txt", **config)[0]
    params:
        atac_samples=lambda wildcards: get_samples(atac_samples, config["contrasts"], wildcards.condition),
        genome=config["genome"],
        jaccard=config["jaccard"],
    threads: 1
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir
        
        ananse binding \
        -A {input.atac} \
        -c {params.atac_samples} \
        -g {params.genome} \
        -p {input.pfm} \
        --pfmscorefile {input.pfmscorefile} \
        --jaccard-cutoff {params.jaccard} \
        -n {threads} \
        -o $outdir \
        > {log} 2>&1
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
        genes=config["rna_tpm"],
    output:
        expand("{result_dir}/network/{{condition}}.tsv",**config),
    log:
        expand("{result_dir}/network/log_{{condition}}.txt",**config),
    benchmark:
        expand("{result_dir}/benchmarks/network_{{condition}}.txt",**config)[0]
    params:
        rna_samples=lambda wildcards: get_samples(rna_samples, config["contrasts"], wildcards.condition),
        genome=config["genome"],
    threads: 1
    resources:
        mem_gb=24
    conda: "../envs/ananse.yaml"
    shell:
        """
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        ananse network \
        {input.binding} \
        -e {input.genes} \
        -c {params.rna_samples} \
        -g {params.genome} \
        -o {output} \
        --full-output \
        -n {threads} \
        > {log} 2>&1
        """


def get_conditions(wildcards):
    networks = dict()
    batch, column, target, source = parse_contrast(wildcards.contrast, rna_samples, check=True)
    networks["target"] = f"{config['result_dir']}/network/{target}.tsv"
    networks["source"] = f"{config['result_dir']}/network/{source}.tsv"
    return networks

rule influence:
    """
    Find the most influential TFs between to conditions
    - as defined in the config.yaml, under contrasts.
    """
    input:
        unpack(get_conditions),
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
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

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
        outdir=$(dirname {output})

        # for the log
        mkdir -p $outdir

        ananse plot \
        {input.inf} \
        -d {input.diff_inf} \
        -t {params.type} \
        --full-output \
        -o $outdir \
        > {log} 2>&1
        """
