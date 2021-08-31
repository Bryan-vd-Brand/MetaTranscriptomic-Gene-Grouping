#This ruleset contains the 3rd step of the pipeline, selecting crAssPhage genomes for further analysis.
#This selection is based on the % amount of mapping to the different available genomes


GENOMES_DIR = config['genomes_dir']
COMBINED_REFERENCE_FILE = config['combined_reference_file']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
    
def get_pileups():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        idxstatFile = f"results/4_ReferenceSelection/per_sample/{RunAccession}/pileup_{RunAccession}.txt"
        results.append(idxstatFile)
    return results

ruleorder: STAR_align_paired > STAR_align_single


rule STAR_index_generation:
    input:
        srcdir(COMBINED_REFERENCE_FILE)
    output:
        "finished_index.touch"
    shell:
        """
        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {GENOMES_DIR} --genomeFastaFiles {input} --genomeSAindexNbases 7
        touch finished_index.touch
        """

rule STAR_align_single:
    input:
        fq = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz",
        indexed = "finished_index.touch"
    output:
        "results/4_ReferenceSelection/per_sample/{sample}/{sample}.idxstat",
        "results/4_ReferenceSelection/per_sample/{sample}/pileup_{sample}.txt"
    shell:
        """
        STAR --runThreadN 20 --alignIntronMax 1 --genomeDir {GENOMES_DIR} --readFilesIn {input.fq} --readFilesCommand gunzip -c --outFileNamePrefix {wildcards.sample} ;
        samtools view -O BAM -b {wildcards.sample}Aligned.out.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/4_ReferenceSelection/per_sample/{wildcards.sample}/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/4_ReferenceSelection/per_sample/{wildcards.sample}/{wildcards.sample}.flagstat ;
        pileup.sh overwrite=true in={wildcards.sample}Aligned.out.sam out=results/4_ReferenceSelection/per_sample/{wildcards.sample}/pileup_{wildcards.sample}.txt ;
        rm -rf {wildcards.sample}.bam ;
        rm -rf {wildcards.sample}_mapped.bam ;
        rm -rf {wildcards.sample}_mapped.bam.bai ;
        rm -rf sorted_{wildcards.sample}.bam ;
        rm -rf {wildcards.sample}_mapped_q3.bam ;
        rm -rf {wildcards.sample}_mapped_q3.bam.bai ;
        rm -rf {wildcards.sample}Aligned.out.sam ;
        rm -rf {wildcards.sample}.bam.bai ;
        rm -df {wildcards.sample}_STARtmp ;
        mv {wildcards.sample}*.out ./results/4_ReferenceSelection/per_sample/{wildcards.sample}/ ;
        """

rule STAR_align_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz",
        indexed = "finished_index.touch"
    output:
        "results/4_ReferenceSelection/per_sample/{sample}/{sample}.idxstat",
        "results/4_ReferenceSelection/per_sample/{sample}/pileup_{sample}.txt"
    shell:
        """
        STAR --runThreadN 20 --alignIntronMax 1 --genomeDir {GENOMES_DIR} --readFilesIn {input.r1} {input.r2} --readFilesCommand gunzip -c --outFileNamePrefix {wildcards.sample} ;
        samtools view -O BAM -b {wildcards.sample}Aligned.out.sam > {wildcards.sample}.bam ;
        samtools sort {wildcards.sample}.bam > sorted_{wildcards.sample}.bam ;
        samtools view -b -F 4 sorted_{wildcards.sample}.bam > {wildcards.sample}_mapped.bam ;
        samtools view -O BAM -q 3 {wildcards.sample}_mapped.bam > {wildcards.sample}_mapped_q3.bam ;
        samtools index {wildcards.sample}_mapped_q3.bam ;
        samtools idxstats {wildcards.sample}_mapped_q3.bam > results/4_ReferenceSelection/per_sample/{wildcards.sample}/{wildcards.sample}.idxstat ;
        samtools flagstat -O tsv {wildcards.sample}_mapped_q3.bam > results/4_ReferenceSelection/per_sample/{wildcards.sample}/{wildcards.sample}.flagstat ;
        pileup.sh overwrite=true in={wildcards.sample}Aligned.out.sam out=results/4_ReferenceSelection/per_sample/{wildcards.sample}/pileup_{wildcards.sample}.txt ;
        rm -rf {wildcards.sample}.bam ;
        rm -rf {wildcards.sample}_mapped.bam ;
        rm -rf {wildcards.sample}_mapped.bam.bai ;
        rm -rf sorted_{wildcards.sample}.bam ;
        rm -rf {wildcards.sample}_mapped_q3.bam ;
        rm -rf {wildcards.sample}_mapped_q3.bam.bai ;
        rm -rf {wildcards.sample}Aligned.out.sam ;
        rm -rf {wildcards.sample}.bam.bai ;
        rm -df {wildcards.sample}_STARtmp ;
        mv {wildcards.sample}*.out ./results/4_ReferenceSelection/per_sample/{wildcards.sample}/ ;
        """

rule select_genomes:
    input:
        pileups = get_pileups()
    output:
        out_table  = report("results/4_ReferenceSelection/ReferenceTable.tsv", category="References"),
        out_list = report('results/4_ReferenceSelection/ReferenceList.txt', category="References"),
    params:
        script = srcdir("../Scripts/gather_references.py")
    shell:
        "python {params.script} -i {input.pileups} "
        
rule generate_hcov_heatmap:
    input:
        "results/4_ReferenceSelection/ReferenceTable.tsv"
    params:
        script = srcdir("../Scripts/generateHcovHeatmap.R")
    output:
        "results/4_ReferenceSelection/finished_plotgeneration.touch"
    shell:
        """
        Rscript {params.script} \
        ./results/4_ReferenceSelection/per_sample ;
        touch ./results/4_ReferenceSelection/finished_plotgeneration.touch ;
        """
