#This snakefile contains a rule that takes pileup results generated in 5_ and creates 2 histogram showing the the hcov&stdDev distribution for the selected genomes    

rule Generate_Horizontal_Coverage_Plot:
    input:
        "results/5_SingleReferenceGeneExpression/finished_allsamples.touch"
    params:
        script = srcdir("../Scripts/generateHCovPlot.R")
    output:
        "results/6_GenerateHorizontalCoveragePlot/finished_plotgeneration.touch"
    shell:
        """
        Rscript {params.script} \
        ./results/5_SingleReferenceGeneExpression/per_sample ;
        touch ./results/6_GenerateHorizontalCoveragePlot/finished_plotgeneration.touch ;
        """