#This rule file performs fastqc/multiqc on the original data set.

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
	
def get_all_samples():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        results.append(RunAccession)
    return results

ruleorder: fastqc_paired > fastqc_single

rule fastqc_single:
    input:
        f"{DATA_DIR}/{{sample}}/{{sample}}.fastq.gz"
    output:
        multiext("results/0_fastqc/{sample}/{sample}_fastqc", ".zip", ".html"),
        outdir = directory("results/0_fastqc/{sample}")
    threads: 2
    log:
        stdout = "results/0_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "results/0_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """

rule fastqc_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/{{sample}}_2.fastq.gz"
    output:
        expand("results/0_fastqc/{{sample}}/{{sample}}_{read}_fastqc.zip", read=["1", "2"]),
        expand("results/0_fastqc/{{sample}}/{{sample}}_{read}_fastqc.html", read=["1", "2"]),
        outdir = directory("results/0_fastqc/{sample}")
    threads: 2
    log:
        stdout = "results/0_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "results/0_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """

rule multiqc:
    input:
        expand("results/0_fastqc/{sample}", sample=get_all_samples())
    output:
        "results/0_fastqc/multiqc_report.html"
    params:
        outdir = "results/0_fastqc"
    shell:
        """
        multiqc -f {input} -o {params.outdir}
        """
