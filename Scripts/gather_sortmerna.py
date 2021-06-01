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
    requiredArgs.add_argument(
        "-d",
        "--input_dna_files",
        dest = "input_dna_files",
        nargs="+",
        required=True,
        help=""
    )
    requiredArgs.add_argument(
        "-m",
        "--mock_samples",
        dest = "mock_samples",
        nargs="+",
        required=True,
        help=""
    )

    return parser.parse_args()

def gather_sortmerna(sortmerna_files, mock_samples):
    samples = list()
    percentages = list()
    transcr = list()

    for file in sortmerna_files:
        sample = os.path.basename(file).split("_sortmerna")[0]
        samples.append(sample)

        if sample not in mock_samples:
            transcr.append("yes")
        else:
            transcr.append("no")

        lines = [line.strip() for line in open(file).readlines()]
        for line in lines:
            if "Total reads passing E-value threshold" in line:
                percentage = float(line.split("(")[1].split(")")[0])
                percentages.append(percentage)

    df = pd.DataFrame(list(zip(samples, percentages, transcr)),
               columns =['samples', 'rRNA percent', 'transcr'])
    df = df.sort_values("rRNA percent")

    palette = {"yes": "blue", "no":"red"}
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    fig.patch.set_facecolor('white')
    sns.barplot(x="samples", y="rRNA percent", data=df, hue="transcr",dodge=False, palette=palette)
    plt.xticks(rotation=60)

    for (type, ticklbl) in zip(df['transcr'], ax.xaxis.get_ticklabels()):
        ticklbl.set_color('blue' if type == 'yes' else 'red')

    fig.savefig("results/0_sortmerna/rRNA_percentages.png")

    df.to_csv("results/0_sortmerna/rRNA_percentages.txt", sep="\t", index=False)

def main():

    args = parse_args()

    all_files = args.input_files + args.input_dna_files

    gather_sortmerna(all_files, args.mock_samples)

if __name__ == "__main__":
    main()
