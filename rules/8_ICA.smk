#This snakemake rule file contains rules associated with starting the Independent Component Analysis from the expression counts of all the genes. 

VARIANCE = config['pca_variance_required']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split

rule generate_distribution_graph:
    input:
        "results/8_ICA/finished_ICA.touch"
    params:
        script = srcdir("../Scripts/generateDistribution.R")
    output:
        "results/8_ICA/generate_distribution_graph.touch"
    shell:
        """
        Rscript {params.script} \
        ./results/8_ICA ;
        touch ./results/8_ICA/generate_distribution_graph.touch ;
        """

rule Run_ICA:
    input:
        "FeatureCount_table.tsv"
    params:
        script = srcdir("../Scripts/runICA.py")
    output:
        "results/8_ICA/finished_ICA.touch"
    shell:
        """
        python {params.script} -ft {input} -PCAvar {VARIANCE} ;
        touch ./results/8_ICA/finished_ICA.touch ;
        """
