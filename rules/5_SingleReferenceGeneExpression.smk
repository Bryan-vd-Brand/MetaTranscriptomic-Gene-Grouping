#This snakemake file contains the rules associated with generating featurecount results for sample-genome combinations

GENE_ANNOTATION_FILE = config.get("gene_annotation_file")
GENOME_DIR = config['seperated_genomes_dir']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
    
def get_references():
    allReferences = []
    if not Path('./results/4_ReferenceSelection/ReferenceList.txt').exists():
        return []
    with open('results/4_ReferenceSelection/ReferenceList.txt', 'r') as referenceList:
        for ref in referenceList:
            allReferences.append(ref.rstrip().split('\t')[0])
    return allReferences
    
def get_all_samples():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        results.append(RunAccession)
    return results

checkpoint check_referenceList:
    input:
          "results/4_ReferenceSelection/ReferenceList.txt"
    output:
        touch(".check_ref.touch") #fakeoutput?

class CheckPoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, w):
        global checkpoints
        checkpoints.check_referenceList.get(**w)
        references = get_references()
        pattern = expand(F"{GENOME_DIR}/{{ref}}.fasta", ref = references)
        return pattern

ruleorder: run_genomes_paired > run_genomes_single

rule index_genomes:
    input:
        refsCount = "results/4_ReferenceSelection/ReferenceList.txt",
        refs = expand(F"{GENOME_DIR}/{{ref}}.fasta", ref = get_references()),
    output:
        "results/5_SingleReferenceGeneExpression/finished_Index.touch"
    params:
        script = srcdir("../Scripts/make_topGenome_indexes.py")
    shell:
        "python {params.script} -rl {input.refsCount} -g {input.refs} -gd {GENOME_DIR}"

rule run_genomes_single:
    input:
        sample = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz",
        refs = CheckPoint_MakePattern(""),
        refsCount = "results/4_ReferenceSelection/ReferenceList.txt",
        refsTable = "results/4_ReferenceSelection/ReferenceTable.tsv",
        finishedIndex = "results/5_SingleReferenceGeneExpression/finished_Index.touch"
    output:
        "results/5_SingleReferenceGeneExpression/per_sample/finished_{sample}.touch"
    params:
        script = srcdir("../Scripts/run_genomes.py")
    shell:
        "python {params.script} -ra {wildcards.sample} -s {input.sample} -g {input.refs} -a {GENE_ANNOTATION_FILE} -rl {input.refsCount} -rt {input.refsTable} -gd {GENOME_DIR}"

rule run_genomes_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz",
        refs = CheckPoint_MakePattern(""),
        refsCount = "results/4_ReferenceSelection/ReferenceList.txt",
        refsTable = "results/4_ReferenceSelection/ReferenceTable.tsv",
        finishedIndex = "results/5_SingleReferenceGeneExpression/finished_Index.touch"
    output:
        "results/5_SingleReferenceGeneExpression/per_sample/finished_{sample}.touch"
    params:
        script = srcdir("../Scripts/run_genomes.py")
    shell:
        "python {params.script} -ra {wildcards.sample} -s {input.r1} {input.r2} -g {input.refs} -a {GENE_ANNOTATION_FILE} -rl {input.refsCount} -rt {input.refsTable} -gd {GENOME_DIR}"
        
rule finished_genomes:
    input:
        expand("results/5_SingleReferenceGeneExpression/per_sample/finished_{sample}.touch", sample = get_all_samples())
    output:
        "results/5_SingleReferenceGeneExpression/finished_allsamples.touch"
    shell:
        """
        touch ./results/5_SingleReferenceGeneExpression/finished_allsamples.touch ;
        """
