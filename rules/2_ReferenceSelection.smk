#This ruleset contains the 3rd step of the pipeline, selecting crAssPhage genomes for further analysis.
#This selection is based on the % amount of mapping to the diffferent available genomes
#Genomes are sorted by # [unique??] hits for the sample, selecting the top genomes until atleast 50% of the mapped reads are represented (N50 like)
#TODO: Trim adapters
#TODO: Align fastq against genomes, samtools steps, create best genome list, select best, store best

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
    
def get_idxstats():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        idxstatFile = f"results/2_ReferenceSelection/per_sample/{RunAccession}.idxstat"
        results.append(idxstatFile)
    return results

ruleorder: fastp_pe > fastp_se

rule fastp_se:
    input:
        sample=[f"{DATA_DIR}/{{sample}}/{{sample}}.fastq.gz"]
    output:
        trimmed=f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz",
        html="results/2_ReferenceSelection/logs/{sample}.html",
        json="results/2_ReferenceSelection/logs/{sample}.json"
    log:
        "results/2_ReferenceSelection/logs_fastp/{sample}.log"
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
        html="results/2_ReferenceSelection/logs/{sample}.html",
        json="results/2_ReferenceSelection/logs/{sample}.json"
    log:
        "results/2_ReferenceSelection/logs_fastp/{sample}.log"
    params:
        extra=""
    threads: 2
    wrapper:
        "v0.75.0/bio/fastp"

ruleorder: bwa_mem_align_paired > bwa_mem_align_single

rule bwa_mem_align_single:
    input:
        f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz"
    output:
        "results/2_ReferenceSelection/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/2_ReferenceSelection/logs/{sample}_RefSelec.stdout",
        stderr = "results/2_ReferenceSelection/logs/{sample}_RefSelec.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input} > {wildcards.sample}.sam ;
        samtools view -O BAM -b {wildcards.sample}.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{wildcards.sample}.flagstat ;
        pileup.sh in={wildcards.sample}.sam out=results/2_ReferenceSelection/per_sample/pileup_{wildcards.sample}.txt ;
        """

rule bwa_mem_align_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz"
    output:
        "results/2_ReferenceSelection/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/2_ReferenceSelection/logs/{sample}_RefSelec.stdout",
        stderr = "results/2_ReferenceSelection/logs/{sample}_RefSelec.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input.r1} {input.r2} > {wildcards.sample}.sam ;
        samtools view -O BAM -b {wildcards.sample}.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{wildcards.sample}.flagstat ;
        pileup.sh in={wildcards.sample}.sam out=results/2_ReferenceSelection/per_sample/pileup_{wildcards.sample}.txt ;
        """

rule select_genomes:
    input:
        idxfiles = get_idxstats()
    output:
        out_table  = report("results/2_ReferenceSelection/ReferenceTable.tsv", category="References"),
    params:
        script = srcdir("../Scripts/gather_references.py")
    shell:
        "python {params.script} -i {input.idxfiles} "
