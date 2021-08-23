#This snakemake rule file contains rules associated with aggregating all the featurecount results for each of the genomes, creates a table for them and generates a heatmap plot for visual inspection

GENOMES_DIR = config['genomes_dir']
COMBINED_REFERENCE_FILE = config['combined_reference_file']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
	
	
	
rule Generate_FeatureCount_Tables:
    input:
        "results/5_SingleReferenceGeneExpression/finished_allsamples.touch"
    params:
        script = srcdir("../Scripts/gather_FeatureCounts.py")
    output:
        "FeatureCount_table.tsv"
    shell:
        "python {params.script} -rf results/5_SingleReferenceGeneExpression/per_sample/ -rl results/4_ReferenceSelection/ReferenceList.txt"

rule Generate_Heatmap:
    input:
	    "FeatureCount_table.tsv"
    params:
        script = srcdir("../Scripts/generateHeatmap.R")
    output:
        "results/7_FeatureCounts/finished_plotgeneration.touch"
    shell:
        """
        Rscript {params.script} \
        FeatureCount_table.tsv ;
        touch ./results/7_FeatureCounts/finished_plotgeneration.touch ;
        """