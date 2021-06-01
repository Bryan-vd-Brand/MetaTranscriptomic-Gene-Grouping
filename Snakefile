configfile: "config/config.yaml"
report: "report/workflow.rst"
include: "rules/1_sortmerna.smk"
include: "rules/0_fastqc.smk"


# both methods return the same information. Notice how for the paired end samples,
# there is only one string containing the two files. Because of that we have to
# split by "," later to get the two files
print(config["samples"])
print(config.get("samples"))

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split


print(config.get('input_dir'))
DATA_DIR = config.get('input_dir')

DNA_SAMPLES = config.get("metagenomic_DNA_samples")
SORTMERNA_DB_FILES = config.get("sortmerna_db_files")

# Create mock empty input files for DNA_samples
os.makedirs("resources/DNA_samples", exist_ok = True)
for dna_sample in DNA_SAMPLES:
    os.system(f"touch resources/DNA_samples/{dna_sample}.mock")

rule all:
    input:
        "results/0_fastqc/multiqc_report.html",
        "results/1_sortmerna/{sample}_sortmerna.txt"
        #"results/1_sortmerna/rRNA_percentages.png",
        #"results/1_sortmerna/rRNA_percentages.txt",