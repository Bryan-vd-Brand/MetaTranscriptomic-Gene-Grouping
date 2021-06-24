import os
import pandas as pd
import argparse

AnnotationFile = "/home/bryan/Documents/Git/MetaTranscriptomic-Gene-Grouping/resources/genes/seperatedCrAssPhageGeneAnnotation.table"

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-s",
        "--sample_files",
        dest = "sample_files",
        nargs="+",
        required=True,
        help="List of fastq sample names"
    )
    
    requiredArgs.add_argument(
        "-r",
        "--genome_files",
        dest = "genome_files",
        nargs="+",
        required=True,
        help="List of Fasta"
    )
    
    return parser.parse_args()
    
def create_SingleRef_Runs(samples, genomes):
    print("Test")
    print(samples)
    print(genomes)
    for genome in genomes:
        if len(samples) > 1:
            #paired
            bwa_mem_paired_fc(samples, genome) 
        else:
            #unpaired
            bwa_mem_single_fc(samples, genome)
    

#f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz"
#TODO: Fix export file names/structure!
#TODO: Add generation of index files for the genome reference (and deletion??)
def bwa_mem_single_fc(sample, genome):
    RunAccession = sample[0].split("_")[1].split(".")[0]
    print(F"Running {genome} and {RunAccession}")
    os.system(F"bwa mem -t 5 {genome} {sample} > {RunAccession}.sam")
    os.system(F"samtools view -O BAM -b {RunAccession}.sam > {RunAccession}.bam")
    os.system(F"samtools sort {RunAccession}.bam > sorted_{RunAccession}.bam")
    os.system(F"samtools view -b -F 4 sorted_{RunAccession}.bam > {RunAccession}_mapped.bam")
    os.system(F"samtools index {RunAccession}_mapped.bam")
    os.system(F"samtools idxstats {RunAccession}_mapped.bam > results/3_SingleReferenceGeneExpression/per_sample/{RunAccession}.idxstat")
    os.system(F"samtools flagstat -O tsv {RunAccession}_mapped.bam > results/3_SingleReferenceGeneExpression/per_sample/{RunAccession}.flagstat")
    os.system(F"pileup.sh in={RunAccession}.sam out=results/3_SingleReferenceGeneExpression/per_sample/pileup_{RunAccession}.txt")
    os.system(F"featureCounts -M -O -G {genome} -F SAF -a {AnnotationFile} -o results/3_SingleReferenceGeneExpression/per_sample/{RunAccession} {RunAccession}_mapped.bam")
    os.system(F"rm -rf *.bam")
    os.system(F"rm -rf *.sam")
    os.system(F"rm -rf *.bam.bai")

def bwa_mem_paired_fc(sample, genome):
    RunAccession = sample[0].split("_")[1]
    print(F"Running {genome} and {RunAccession}")
    os.system(F"bwa mem -t 5 {genome} {sample[0]} {sample[1]} > {RunAccession}.sam")
    os.system(F"samtools view -O BAM -b {RunAccession}.sam > {RunAccession}.bam")
    os.system(F"samtools sort {RunAccession}.bam > sorted_{RunAccession}.bam")
    os.system(F"samtools view -b -F 4 sorted_{RunAccession}.bam > {RunAccession}_mapped.bam")
    os.system(F"samtools index {RunAccession}_mapped.bam")
    os.system(F"samtools idxstats {RunAccession}_mapped.bam > results/3_SingleReferenceGeneExpression/per_sample/{RunAccession}.idxstat")
    os.system(F"samtools flagstat -O tsv {RunAccession}_mapped.bam > results/3_SingleReferenceGeneExpression/per_sample/{RunAccession}.flagstat")
    os.system(F"pileup.sh in={sample}.sam out=results/3_SingleReferenceGeneExpression/per_sample/pileup_{RunAccession}.txt")
    os.system(F"featureCounts -M -O -G {genome} -F SAF -a {AnnotationFile} -o results/3_SingleReferenceGeneExpression/per_sample/{RunAccession} {RunAccession}_mapped.bam")
    os.system(F"rm -rf *.bam")
    os.system(F"rm -rf *.sam")
    os.system(F"rm -rf *.bam.bai")

def main():

    args = parse_args()
    create_SingleRef_Runs(args.sample_files, args.genome_files)

if __name__ == "__main__":
    main()