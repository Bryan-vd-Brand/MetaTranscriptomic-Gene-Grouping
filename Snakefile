configfile: "config/config.yaml"
report: "report/workflow.rst"
include: "rules/0_fastqc.smk"
include: "rules/1_readTrimming.smk"
include: "rules/2_fastqc.smk"
include: "rules/3_sortmerna.smk"
include: "rules/4_ReferenceSelection.smk"
include: "rules/5_SingleReferenceGeneExpression.smk"
include: "rules/6_GenerateHorizontalCoveragePlot.smk"
include: "rules/7_FeatureCounts.smk"
include: "rules/8_ICA.smk"

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

rule all:
    input:
        "results/0_fastqc/multiqc_report.html",
        "results/2_fastqc/multiqc_report.html",
        "results/3_sortmerna/rRNA_percentages.png",
        "results/3_sortmerna/rRNA_percentages.txt",
        "results/4_ReferenceSelection/ReferenceTable.tsv",
        "results/4_ReferenceSelection/finished_plotgeneration.touch",
        "results/5_SingleReferenceGeneExpression/finished_allsamples.touch",
        "results/6_GenerateHorizontalCoveragePlot/finished_plotgeneration.touch",
        "results/7_FeatureCounts/finished_plotgeneration.touch",
        "results/8_ICA/generate_distribution_graph.touch",
        "results/8_ICA/finished_ICA.touch",
        "results/8_ICA/generated_GeneModules.touch"
