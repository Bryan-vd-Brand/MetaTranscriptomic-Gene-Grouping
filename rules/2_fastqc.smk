#This rule file performs fastqc/multiqc after trimming

DATA_DIR = config['input_dir']

"""
This function is "special" because of the "wildcards" argument. They are called
'input functions' and are briefly describe here https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions .
In short, you can define functions that work with the wildcards as the rules do.
When calling the function in the rule, there is no need to specify the wildcard
as an argument. THEY HAVE A SINGLE ARGUMENT (wilcards), you can not add more
"""
def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split

ruleorder: trimmed_fastqc_paired > trimmed_fastqc_single

rule trimmed_fastqc_single:
    input:
        f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz"
    output:
        multiext("results/2_fastqc/{sample}/trimmed_{sample}_fastqc", ".zip", ".html"),
        outdir = directory("results/2_fastqc/{sample}")
    threads: 2
    log:
        stdout = "results/2_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "results/2_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """

rule trimmed_fastqc_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz"
    output:
        expand("results/2_fastqc/{{sample}}/trimmed_{{sample}}_{read}_fastqc.zip", read=["1", "2"]),
        expand("results/2_fastqc/{{sample}}/trimmed_{{sample}}_{read}_fastqc.html", read=["1", "2"]),
        outdir = directory("results/2_fastqc/{sample}")
    threads: 2
    log:
        stdout = "results/2_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "results/2_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """

rule trimmed_multiqc:
    input:
        expand("results/2_fastqc/{sample}", sample=config.get("samples").keys())
    output:
        "results/2_fastqc/multiqc_report.html"
    params:
        outdir = "results/2_fastqc"
    shell:
        """
        multiqc -f {input} -o {params.outdir}
        """
