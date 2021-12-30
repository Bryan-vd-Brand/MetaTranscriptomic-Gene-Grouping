#This snakemake rule file contains rules associated with starting the Independent Component Analysis from the expression counts of all the genes. 

SAF_ANNOT = config['saf_gene_annotation_file']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split
        
rule perm_statistics:
    input:
        GroupedGeneModules = "results/8_ICA/GroupedGeneModule.tsv"
    params:
        script = srcdir("../Scripts/PermutationStatistics.py")
    output:
        "results/9_PermutationStatistics/GeneModulePval.tsv"
    shell:
        """
        python {params.script} -gaf {SAF_ANNOT} -ggm {input.GroupedGeneModules}
        touch ./results/9_PermutationStatistics/permutation_statistics.touch ;
        """

rule visualize_gene_modules:
    input:
        GeneModules = "results/8_ICA/GeneModules.tsv",
        FC_Table = "FeatureCount_table.tsv",
        PvalGeneModules = "results/9_PermutationStatistics/GeneModulePval.tsv"
    params:
        script = srcdir("../Scripts/visualize_gene_modules.R")
    output:
        "results/9_PermutationStatistics/visualize_gene_modules.touch"
    shell:
        """
        Rscript {params.script} \
        {input.FC_Table} \
        {input.GeneModules} \
        {input.PvalGeneModules} \
        {SAF_ANNOT} ;
        touch ./results/9_PermutationStatistics/visualize_gene_modules.touch ;
        """

rule visualize_ICA_P:
    input:
        PvalGeneModules = "results/9_PermutationStatistics/GeneModulePval.tsv"
    params:
        script = srcdir("../Scripts/visualize_ICA_P.R")
    output:
        "results/9_PermutationStatistics/visualize_ICA_P.touch"
    shell:
        """
        Rscript {params.script} \
        {input.PvalGeneModules} ;
        touch ./results/9_PermutationStatistics/visualize_ICA_P.touch ;
        """