import glob, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-i",
        "--input_files",
        dest = "input_files",
        nargs="+",
        required=True,
        help=""
    )
    return parser.parse_args()

def gather_sortmerna(sortmerna_files):
    samples = list()
    percentages = list()

    for file in sortmerna_files:
        sample = os.path.basename(file).split("_sortmerna")[0]
        samples.append(sample)
        lines = [line.strip() for line in open(file).readlines()]
        for line in lines:
            if "Total reads passing E-value threshold" in line:
                percentage = float(line.split("(")[1].split(")")[0])
                percentages.append(percentage)

    df = pd.DataFrame(list(zip(samples, percentages)),
               columns =['samples', 'rRNA percent'])
    df = df.sort_values("rRNA percent")

    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    fig.patch.set_facecolor('white')
    sns.barplot(x="samples", y="rRNA percent", data=df,dodge=False)
    plt.xticks(rotation=60)
    
    fig.savefig("results/3_sortmerna/rRNA_percentages.png")
    df.to_csv("results/3_sortmerna/rRNA_percentages.txt", sep="\t", index=False)

def main():

    args = parse_args()
    gather_sortmerna(args.input_files)

if __name__ == "__main__":
    main()
