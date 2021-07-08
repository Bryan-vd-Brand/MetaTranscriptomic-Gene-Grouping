import os
import pandas as pd
import argparse




def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")
    
    requiredArgs.add_argument(
        "-ra",
        "--run_accession",
        dest = "run_accession",
        nargs="+",
        required=True,
        help="RunAccession Identifier of Sample"
    )

    requiredArgs.add_argument(
        "-s",
        "--sample_files",
        dest = "sample_files",
        nargs="+",
        required=True,
        help="List of fastq sample names"
    )
    
    requiredArgs.add_argument(
        "-g",
        "--genome_files",
        dest = "genome_files",
        nargs="+",
        required=True,
        help="List of Fasta"
    )
    
    requiredArgs.add_argument(
        "-a",
        "--AnnotationFileGenes",
        dest = "annotation_file",
        nargs="+",
        required=True,
        help=".table file containing genes in SAF"
    )
    
    requiredArgs.add_argument(
        "-rl",
        "--SelectedReferences",
        dest = "reference_list",
        nargs="+",
        required=True,
        help=".tsv file containing references selected for in previous rule"
    )
    
    requiredArgs.add_argument(
        "-rt",
        "--ReferenceTable",
        dest = "reference_table",
        nargs="+",
        required=True,
        help=".tsv file containing present references for each sample"
    )
    
    requiredArgs.add_argument(
        "-gd",
        "--GenomesDirectory",
        dest = "genomes_directory",
        nargs="+",
        required=True,
        help="Full path to a directory containing .fasta files for each of the genomes used"
    )
    
    
    
    return parser.parse_args()
    


#TODO: Fix export file names/structure!
#TODO: Add generation of index files for the genome reference (and deletion??)
def STAR_single_fc(sample, sampleName, referenceDir, genomeName, annotation):
    print(F"Running {genomeName} and {sampleName}")
    os.system(F"STAR --runThreadN 20 --alignIntronMax 1 --genomeDir {referenceDir} --readFilesIn {sample[0]} --readFilesCommand gunzip -c --outFileNamePrefix {sampleName}")
    os.system(F"samtools view -O BAM -b {sampleName}Aligned.out.sam > {sampleName}.bam")
    os.system(F"samtools sort {sampleName}.bam > sorted_{sampleName}.bam")
    os.system(F"samtools view -b -F 4 sorted_{sampleName}.bam > {sampleName}_mapped.bam")
    os.system(F"samtools index {sampleName}_mapped.bam")
    os.system(F"samtools idxstats {sampleName}_mapped.bam > results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}.idxstat")
    os.system(F"samtools flagstat -O tsv {sampleName}_mapped.bam > results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}.flagstat")
    os.system(F"pileup.sh in={sampleName}Aligned.out.sam out=results/5_SingleReferenceGeneExpression/per_sample/pileup_{genomeName}_{sampleName}.txt")
    os.system(F"featureCounts -M -O -F SAF -a {annotation[0]} -o results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}_featurecount {sampleName}_mapped.bam")
    os.system(F"rm -rf {sampleName}.bam")
    os.system(F"rm -rf sorted_{sampleName}.bam")
    os.system(F"rm -rf {sampleName}_mapped.bam")
    os.system(F"rm -rf {sampleName}_mapped.bam.bai")
    os.system(F"rm -rf {sampleName}Aligned.out.sam")
    os.system(F"rm -rf {sampleName}.bam.bai")
    os.system(F"rm -df {sampleName}_STARtmp")
    os.system(F"mv {sampleName}*.out ./results/5_SingleReferenceGeneExpression/per_sample/")

def STAR_paired_fc(sample, sampleName, referenceDir, genomeName, annotation):
    print(F"Running {genomeName} and {sampleName} from {sample}")
    os.system(F"STAR --runThreadN 20 --alignIntronMax 1 --genomeDir {referenceDir} --readFilesIn {sample[0]} {sample[1]} --readFilesCommand gunzip -c --outFileNamePrefix {sampleName}")
    os.system(F"samtools view -O BAM -b {sampleName}Aligned.out.sam > {sampleName}.bam")
    os.system(F"samtools sort {sampleName}.bam > sorted_{sampleName}.bam")
    os.system(F"samtools view -b -F 4 sorted_{sampleName}.bam > {sampleName}_mapped.bam")
    os.system(F"samtools index {sampleName}_mapped.bam")
    os.system(F"samtools idxstats {sampleName}_mapped.bam > results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}.idxstat")
    os.system(F"samtools flagstat -O tsv {sampleName}_mapped.bam > results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}.flagstat")
    os.system(F"pileup.sh in={sampleName}Aligned.out.sam out=results/5_SingleReferenceGeneExpression/per_sample/pileup_{genomeName}_{sampleName}.txt")
    os.system(F"featureCounts -M -O -F SAF -a {annotation[0]} -o results/5_SingleReferenceGeneExpression/per_sample/{genomeName}_{sampleName}_featurecount {sampleName}_mapped.bam")
    os.system(F"rm -rf {sampleName}.bam")
    os.system(F"rm -rf sorted_{sampleName}.bam")
    os.system(F"rm -rf {sampleName}_mapped.bam")
    os.system(F"rm -rf {sampleName}_mapped.bam.bai")
    os.system(F"rm -rf {sampleName}Aligned.out.sam")
    os.system(F"rm -rf {sampleName}.bam.bai")
    os.system(F"rm -df {sampleName}_STARtmp")
    os.system(F"mv {sampleName}*.out ./results/5_SingleReferenceGeneExpression/per_sample/")

#select set of prefered references and do index for them ; run over all samples, if ref in n50 selection do analysis ; create refname_samplename_featurecount result
def main():

    args = parse_args()
    
    references = pd.read_table(args.reference_list[0] , sep='\t', header = None)
    references = references.sort_values(by=[1], ascending=False)
    #Set of best 5 references across all samples
    SelectedReferences = []
    for x in range(4):
        SelectedReferences.append(references[0][x])
          
    #Check for the set of SelectedReferences if they are present in the current sample and if so run it
    referencesTable = pd.read_table(args.reference_table[0], sep = '\t', header = None)
    for refName in SelectedReferences:
        if refName in referencesTable.loc[referencesTable[0] == F'{args.run_accession[0]}', 2].tolist()[0]:
            print(F"FOUND {refName} for {args.run_accession[0]}")
            #Create Index for selected Reference genome
            for referenceFile in args.genome_files:
                if referenceFile.find(refName) > 0:
                    newGenomeDirectory = F"{args.genomes_directory[0]}/{refName}"
                    #if the directory does not exist the index for this genome has not been generated before. 
                    if not os.path.isdir(newGenomeDirectory):
                        os.system(F"mkdir {newGenomeDirectory}")
                        os.system(F"STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {newGenomeDirectory} --genomeFastaFiles {referenceFile} --genomeSAindexNbases 7")
                    if len(args.sample_files) > 1:
                        #paired
                        STAR_paired_fc(args.sample_files, args.run_accession[0], newGenomeDirectory, refName, args.annotation_file) 
                    else:
                        #unpaired
                        STAR_single_fc(args.sample_files, args.run_accession[0], newGenomeDirectory, refName, args.annotation_file)
                        
    #finished for this sample create touch file 
    os.system(F"touch ./results/5_SingleReferenceGeneExpression/per_sample/finished_{args.run_accession[0]}.touch")

if __name__ == "__main__":
    main()