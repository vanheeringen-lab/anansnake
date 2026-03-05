from snakemake.io import expand


# sample: {"peaks": "/path/to/narrowpeak", "bam": "/path/to/bam"}
narrowpeak_files = {}  # TODO: populate
GENOME = None  # TODO


def get_peakfile(wildcards):
    return narrowpeak_files[wildcards.sample]["peaks"]


rule narrowpeak_summit:
    """
    Convert a narrowpeak file to a "macs2 summits" file.
    """
    input:
        get_peakfile
    # noinspection SmkRedundantComma
    output:
        expand("{peak_dir}/{{sample}}_summits.bed", **config),
    log:
        expand("{log_dir}/narrowpeak_summit_{{sample}}.log", **config),
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{ print $1,$2+$10,$2+$10+1,$4,$9; }}' {input} > {output} 2> {log}
        """


def get_summitfiles(wildcards):  # noqa
    return expand(
        [
            f"{{peak_dir}}/{sample}_summits.bed"
            for sample in narrowpeak_files.keys()
        ],
        **config,
    )


rule combine_peaks:
    """
    Uses gimmemotifs' combine_peaks to "combine" peaks. This finds all peaks close
    together and takes the most significant one as the true peak.
    """
    input:
        summitfiles=get_summitfiles,
        sizes=f"{GENOME}.sizes",  # noqa
        genome=GENOME,  # noqa
    output:
        expand("{peak_dir}/combined_summits.bed", **config),
    log:
        expand("{log_dir}/combine_peaks.log", **config),
    conda:
        "../envs/gimme.yaml"
    params:
        windowsize=2 * config.get("peak_windowsize", 100),
    shell:
        """
        combine_peaks --genome {input.genome} --window {params.windowsize} \
        {input.summitfiles} > {output} 2> {log}
        """


rule bedtools_slop:
    """
    After combine_peaks we end up with just a bed file of summits. We extend all peaks
    to the same width, for a fair comparison between peaks.
    """
    input:
        bedfile=rules.combine_peaks.output,  # noqa
        sizes=f"{GENOME}.sizes",  # noqa
    output:
        expand("{peak_dir}/combined_peaks.bed", **config),
    log:
        expand("{log_dir}/bedtools_slop.log", **config),
    conda:
        "../envs/bedtools.yaml"
    params:
        slop=config.get("slop", 100),
    shell:
        """
        bedtools slop -i {input.bedfile} -g {input.sizes} -b {params.slop} | uniq > {output} 2> {log}
        """


rule samtools_index:
    """
    Create an bam index. 
    Assumes the input bam file is sorted.
    """
    input:
        "{filepath}.bam",
    output:
        "{filepath}.bam.bai",
    params:
        config.get("samtools_index", ""),
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index {params} {input} {output}
        """


def get_bams(wildcards):  # noqa
    return [narrowpeak_files[sample]["bam"] for sample in narrowpeak_files.keys()]


def get_bais(wildcards):
    return [f"{f}.bai" for f in get_bams(wildcards)]


# def get_coverage_table_replicates(file_ext):
#     narrowpeak_files[wildcards.sample]["bam"]
#     def wrapped(wildcards):
#         return expand(
#             [
#                 f"{narrowpeak_files[wildcards.sample]["bam"]}.{file_ext}"
#                 for replicate in narrowpeak_files.keys()
#             ],
#             **config,
#         )
#
#     return wrapped
#
#
# def get_names(wildcards):
#     names = [""]
#     for rep in treps[treps["assembly"] == ORI_ASSEMBLIES[wildcards.assembly]].index:
#         names.append(rep_to_descriptive(rep, brep=False))
#     names = "\t".join(names)
#     return names


rule coverage_table:
    """
    Use gimmemotif's coverage_table to generate a count table with the nr of reads
    under each peak per sample.
    """
    input:
        peaks=rules.bedtools_slop.output,  # noqa
        replicates=get_bams,
        replicate_bai=get_bais,
    output:
        expand("{peak_dir}/peak_counts.tsv", **config),
    log:
        expand("{log_dir}/coverage_table.log", **config),
    conda:
        "../envs/gimme.yaml"
    params:
        peak_width=2 * config.get("slop", 100),  # same width as the upstream files
    resources:
        mem_gb=3,
    threads: 12  # default of the function
    shell:
        """
        echo "# The number of reads under each peak" > {output} 
        coverage_table {input.peaks} {input.replicates} --window {params.peak_width} --nthreads {threads} \
        2> {log} | grep -vE "^#" 2>> {log} 
        """
        # |
        # awk 'BEGIN {{ FS = "@" }} NR==1{{gsub("{wildcards.assembly}-|.samtools-coordinate","",$0)}}; \
        # {{print $0}}' >> {output}
        #
        # # overwrite sample names with descriptive/replicate names
        # sed-i "2s/.*/{params.names}/" {output}
        # """
