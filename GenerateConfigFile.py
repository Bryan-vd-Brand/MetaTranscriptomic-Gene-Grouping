#This script opens the config.yaml, reads out the containing lines until #START
#Then takes the data directory, gathers all of the different fastq files present in there
#Requires the following structure DataFile/SampleName/SampleName.fastq.gz (or paired ; SampleName_1.fastq.gz)
#  ERR2088987: Data/ERR2088987/ERR2088987.fastq.gz
#  SRR2992971: Data/SRR2992971/SRR2992971_1.fastq,Data/SRR2992971/SRR2992971_2.fastq.gz
import glob
import os
import pandas as pd

ConfigFile = "./config/config.yaml"
DataFolder = ""
def main(): 
    with open(F'{ConfigFile}','rt') as configFile:
        line = configFile.readline()
        while("input_dir:" not in line):
            line = configFile.readline()
        DataFolder = configFile.readline().strip()
    with open("./config/temp_config.yaml",'w') as newConfigFile:
        with open(F'{ConfigFile}','rt') as configFile:
            line = configFile.readline()
            while("#START" not in line):
                newConfigFile.write(line)
                line = configFile.readline()
            newConfigFile.write(line) #START line, add it again
            #generate new entry for all samples in DataFolder
            newConfigFile.write("samples:\n")
            for sampleFolder in os.listdir(DataFolder):
                if os.path.exists(F"{DataFolder}/{sampleFolder}/{sampleFolder}_1.fastq.gz"): #paired
                    newConfigFile.write(F"  {sampleFolder}: {DataFolder}/{sampleFolder}/{sampleFolder}_1.fastq.gz,{DataFolder}/{sampleFolder}/{sampleFolder}_2.fastq.gz\n")
                elif os.path.exists(F"{DataFolder}/{sampleFolder}/{sampleFolder}.fastq.gz"): #unpaired
                    newConfigFile.write(F"  {sampleFolder}: {DataFolder}/{sampleFolder}/{sampleFolder}.fastq.gz\n")
                else:
                    print(F"ERROR failed to create config entry for {sampleFolder}")
    os.system(F"rm -f {ConfigFile}")
    os.rename("./config/temp_config.yaml","./config/config.yaml")
    print(DataFolder)
    print(F"Finished generating new config file out of ^")

if __name__ == "__main__":
	main()