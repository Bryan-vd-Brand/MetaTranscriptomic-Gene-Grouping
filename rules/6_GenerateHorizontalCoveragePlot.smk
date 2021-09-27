#This snakefile contains a rule that takes pileup results generated in 5_ and creates 2 histogram showing the the hcov&stdDev distribution for the selected genomes    

rule Generate_Horizontal_Coverage_Plot:
    input:
        finished = "results/4_ReferenceSelection/finished_plotgeneration.touch",
        ReferenceList = "results/4_ReferenceSelection/ReferenceList.txt"
    params:
        script = srcdir("../Scripts/generateHCovPlot.R")
    output:
        "results/6_GenerateHorizontalCoveragePlot/finished_plotgeneration.touch"
    shell:
        """
        Rscript {params.script} \
        ./results/4_ReferenceSelection/per_sample \
        {input.ReferenceList} ;
        touch ./results/6_GenerateHorizontalCoveragePlot/finished_plotgeneration.touch ;
        """
        
