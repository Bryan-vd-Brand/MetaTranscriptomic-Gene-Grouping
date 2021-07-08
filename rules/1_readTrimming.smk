#Rule file trims fastq files for adapters and bad sequences using fastp.

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
    
ruleorder: fastp_pe > fastp_se

rule fastp_se:
    input:
        sample=[f"{DATA_DIR}/{{sample}}/{{sample}}.fastq.gz"]
    output:
        trimmed=f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz",
        html="results/1_readTrimming/logs/{sample}.html",
        json="results/1_readTrimming/logs/{sample}.json"
    log:
        "results/1_readTrimming/logs_fastp/{sample}.log"
    params:
        extra=""
    threads: 1
    wrapper:
        "v0.75.0/bio/fastp"


rule fastp_pe:
    input:
        sample=[f"{DATA_DIR}/{{sample}}/{{sample}}_1.fastq.gz", f"{DATA_DIR}/{{sample}}/{{sample}}_2.fastq.gz"]
    output:
        trimmed=[f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz", f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz"],
        html="results/1_readTrimming/logs/{sample}.html",
        json="results/1_readTrimming/logs/{sample}.json"
    log:
        "results/1_readTrimming/logs_fastp/{sample}.log"
    params:
        extra=""
    threads: 2
    wrapper:
        "v0.75.0/bio/fastp"