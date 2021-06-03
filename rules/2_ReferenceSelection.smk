#This ruleset contains the 3rd step of the pipeline, selecting crAssPhage genomes for further analysis.
#This selection is based on the % amount of mapping to the diffferent available genomes
#Genomes are sorted by # [unique??] hits for the sample, selecting the top genomes until atleast 50% of the mapped reads are represented (N50 like)
#TODO: Trim adapters
#TODO: Align fastq against genomes, samtools steps, create best genome list, select best, store best

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split

ruleorder: bwa_mem_align_paired > bwa_mem_align_single

rule bwa_mem_align_single:
    input:
        f"{DATA_DIR}/{{sample}}/{{sample}}.fastq.gz"
    output:
        "results/2_ReferenceSelection/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/2_ReferenceSelection/{sample}/{sample}_fastqc.stdout",
        stderr = "results/2_ReferenceSelection/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input} > {sample}.sam ;
        samtools view -h -b -S {sample}.sam > {sample}.bam ;
        samtools view -b -F 4 {sample}.bam > {sample}_mapped.bam ;
        samtools view -q 3 {sample}_mapped.bam > {sample}_mapped_q3.bam ;
        samtools index {sample}_mapped_q3.bam ;
        samtools idxstats {sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{sample}.idxstat ;
        """

rule bwa_mem_align_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/{{sample}}_2.fastq.gz"
    output:
        "results/2_ReferenceSelection/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/2_ReferenceSelection/{sample}/{sample}_fastqc.stdout",
        stderr = "results/2_ReferenceSelection/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input.r1} {input.r2} > {sample}.sam ;
        samtools view -h -b -S {sample}.sam > {sample}.bam ;
        samtools view -b -F 4 {sample}.bam > {sample}_mapped.bam ;
        samtools view -q 3 {sample}_mapped.bam > {sample}_mapped_q3.bam ;
        samtools index {sample}_mapped_q3.bam ;
        samtools idxstats {sample}_mapped_q3.bam > results/2_ReferenceSelection/per_sample/{sample}.idxstat ;
        """

rule select_genomes:
    input:
        idxfile = "{sample}.idxstat"
    output:
        out_table  = report("results/2_ReferenceSelection/ReferenceTable.txt", category="References"),
    params:
        script = TODO GET SCRIPT
