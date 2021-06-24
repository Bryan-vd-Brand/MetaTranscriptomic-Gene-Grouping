#This snakemake file contains the rules associated with 
#creating the Gene Expression Heatmap for a single reference(genome) across the entire dataset

#TODO: align trimmed fastq against single reference, Samtools [view/sort/index], create flagstat (to store mapped #/% & Supplementary #)
#TODO: Take sorted bam run featureCounts across all fastq
#TODO: Take all featurecount results aggregate into one expression table

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
	
def get_idxstats():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        idxstatFile = f"results/3_SingleReferenceGeneExpression/per_sample/{RunAccession}.idxstat"
        results.append(idxstatFile)
    return results
    
def get_references():
    allReferences = []
    with open('results/2_ReferenceSelection/ReferenceList.txt', 'r') as referenceList:
        for ref in referenceList:
            allReferences.append(ref.rstrip())
    return allReferences

checkpoint check_referenceList:
    input:
          "results/2_ReferenceSelection/ReferenceList.txt"
    output:
        touch(".check_ref.touch") #fakeoutput?

class CheckPoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, w):
        global checkpoints
        checkpoints.check_referenceList.get(**w)
        references = get_references()
        pattern = expand("resources/genomes/seperated/{ref}.fasta", ref = references)
        return pattern

ruleorder: run_genomes_paired > run_genomes_single

rule run_genomes_single:
    input:
        sample = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz",
        refs = CheckPoint_MakePattern("")
    output:
        "results/3_SingleReferenceGeneExpression/per_sample/{sample}.idxstat"
    params:
        script = srcdir("../Scripts/run_genomes.py")
    shell:
        "python {params.script} -s {input.sample} -r {input.refs} "

rule run_genomes_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz",
        refs = CheckPoint_MakePattern("")
    output:
        "results/3_SingleReferenceGeneExpression/per_sample/{sample}.idxstat"
    params:
        script = srcdir("../Scripts/run_genomes.py")
    shell:
        "python {params.script} -s [{input.r1},{input.r2}] -r {input.refs} "

rule select_genomes_test:
    input:
        idxfiles = get_idxstats()
    output:
        out_table  = report("results/3_SingleReferenceGeneExpression/ReferenceTable.tsv", category="References"),
        out_list = report('results/3_SingleReferenceGeneExpression/ReferenceList.txt', category="References"),
    params:
        script = srcdir("../Scripts/gather_references.py")
    shell:
        "python {params.script} -i {input.idxfiles} "