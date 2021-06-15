#This snakemake file contains the rules associated with 
#creating the Gene Expression Heatmap for a single reference(genome) across the entire dataset

#TODO: align trimmed fastq against single reference, Samtools [view/sort/index], create flagstat (to store mapped #/% & Supplementary #)
#TODO: Take sorted bam run featureCounts across all fastq
#TODO: Take all featurecount results aggregate into one expression table

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
    
def get_references():
    allReferences = []
    with open('results/2_ReferenceSelection/ReferenceList.txt', 'r') as referenceList:
        for ref in referenceList:
            allReferences.append(ref.rstrip())
    return allReferences

ruleorder: bwa_mem_OneRef_align_paired > bwa_mem_OneRef_align_single

checkpoint check_referenceList:
    input:
          "results/2_ReferenceSelection/ReferenceList.txt"
    output:
        touch(".check_ref.touch") #fakeoutput?

class CheckPoint_MakePattern:
    def __init__(self,pattern):
        self.pattern = pattern

    def __call__(self, w):
        global checkpoints
        checkpoints.check_referenceList.get(**w)
        references = get_references()
        pattern = expand(self.pattern, ref = references)
        print(references)
        print(pattern)
        return pattern

rule gather_reference_genomes:
    input:
        CheckPoint_MakePattern("resources/genomes/seperated/{ref}.fasta")
    output:
        touch("testfile.touch") #testfile

rule bwa_mem_OneRef_align_single:
    input:
        f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz"
    output:
        "results/3_SingleReferenceGeneExpression/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/3_SingleReferenceGeneExpression/logs/{sample}_RefSelec.stdout",
        stderr = "results/3_SingleReferenceGeneExpression/logs/{sample}_RefSelec.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input} > {wildcards.sample}.sam ;
        samtools view -O BAM -b {wildcards.sample}.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/3_SingleReferenceGeneExpression/per_sample/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/3_SingleReferenceGeneExpression/per_sample/{wildcards.sample}.flagstat ;
        pileup.sh in={wildcards.sample}.sam out=results/3_SingleReferenceGeneExpression/per_sample/pileup_{wildcards.sample}.txt ;
        rm -rf *.bam ;
        rm -rf *.sam ;
        rm -rf *.bam.bai ;
        """

rule bwa_mem_OneRef_align_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz"
    output:
        "results/3_SingleReferenceGeneExpression/per_sample/{sample}.idxstat"
    params:
        index = srcdir("../resources/genomes/crAss_genomes.fasta")
    threads: 5
    log:
        stdout = "results/3_SingleReferenceGeneExpression/logs/{sample}_RefSelec.stdout",
        stderr = "results/3_SingleReferenceGeneExpression/logs/{sample}_RefSelec.stderr"
    shell:
        """
        bwa mem -t 5 {params.index} {input.r1} {input.r2} > {wildcards.sample}.sam ;
        samtools view -O BAM -b {wildcards.sample}.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/3_SingleReferenceGeneExpression/per_sample/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/3_SingleReferenceGeneExpression/per_sample/{wildcards.sample}.flagstat ;
        pileup.sh in={wildcards.sample}.sam out=results/3_SingleReferenceGeneExpression/per_sample/pileup_{wildcards.sample}.txt ;
        rm -rf *.bam ;
        rm -rf *.sam ;
        rm -rf *.bam.bai ;
        """