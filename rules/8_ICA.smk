#This snakemake rule file contains rules associated with starting the Independent Component Analysis from the expression counts of all the genes. 

VARIANCE = config['pca_variance_required']
A_ANNOT = config['additional_annotation_file']
SAF_ANNOT = config['saf_gene_annotation_file']

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

rule generate_gene_modules:
    input:
      "results/8_ICA/finished_ICA.touch"
    params:
        script = srcdir("../Scripts/generateGeneModules.py")
    output:
        "results/8_ICA/GeneModules.tsv"
    shell:
        """
        python {params.script} -rf ./results/8_ICA -gaf {A_ANNOT}
        """
        
rule many_ICA:    
    input:
        "FeatureCount_table.tsv"
    params:
        script = srcdir("../Scripts/1000_runICA.py")
    output:
        "results/8_ICA/GroupedGeneModule.tsv"
    shell:
        """
        python {params.script} -ft {input} -PCAvar {VARIANCE} -gaf {A_ANNOT}
        """
        
rule visualize_gene_modules:
    input:
        GeneModules = "results/8_ICA/GeneModules.tsv",
        FC_Table = "FeatureCount_table.tsv",
        GroupedGeneModules = "results/8_ICA/GroupedGeneModule.tsv"
    params:
        script = srcdir("../Scripts/visualize_gene_modules.R")
    output:
        "results/8_ICA/visualize_gene_modules.touch"
    shell:
        """
        Rscript {params.script} \
        {input.FC_Table} \
        {input.GeneModules} \
        {input.GroupedGeneModules} \
        {SAF_ANNOT} ;
        touch ./results/8_ICA/visualize_gene_modules.touch ;
        """