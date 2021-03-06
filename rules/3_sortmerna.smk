#This rule file runs sortmerna on the provided dataset

import glob

SORTMERNA_DB_FILES = config.get("sortmerna_db_files")
DATA_DIR = config['input_dir']

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split

def get_sortmerna_results():
    results = []
    for file in config["samples"]:
        RunAccession = file.split(":")[0]
        sortmernaFile = f"results/3_sortmerna/{RunAccession}_sortmerna.txt"
        results.append(sortmernaFile)
    return results

ruleorder: run_sortmerna_paired > run_sortmerna_single

rule index_sortmerna:
    input:
        expand("resources/sortmerna_db/{db_file}", db_file=config.get("sortmerna_db_files"))
    output:
        expand("resources/sortmerna_db/{db_file}.done", db_file=SORTMERNA_DB_FILES)
    shell:
        """
        for db_file in resources/sortmerna_db/*fasta;
            do
                sortmerna --ref $db_file --index 1 --idx-dir resources/sortmerna_db > $db_file.log;
                touch $db_file.done ;
            done
        """

rule run_sortmerna_single:
    input:
        index = rules.index_sortmerna.output,
        fastq = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}.fastq.gz"
    output:
        "results/3_sortmerna/{sample}_sortmerna.txt",
    params:
        dbfiles =  [f"--ref resources/sortmerna_db/{dbfile}" for dbfile in SORTMERNA_DB_FILES],
        results_log = "results/3_sortmerna/per_sample/{sample}/out/aligned.log",
    log: "results/3_sortmerna/{sample}_sortmerna.log"
    shell:
        "sortmerna "
        "{params.dbfiles} "
        "--reads {input.fastq} "
        "--threads {threads} "
        "--workdir results/3_sortmerna/per_sample/{wildcards.sample} "
        "--idx-dir resources/sortmerna_db &> {log}; "
        "mv {params.results_log} {output} ; "
        "rm -rf results/3_sortmerna/per_sample/{wildcards.sample}"

rule run_sortmerna_paired:
    input:
        index = rules.index_sortmerna.output,
        r1 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_1.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/trimmed_{{sample}}_2.fastq.gz"
    output:
        "results/3_sortmerna/{sample}_sortmerna.txt",
    params:
        dbfiles =  [f"--ref resources/sortmerna_db/{dbfile}" for dbfile in SORTMERNA_DB_FILES],
        results_log = "results/3_sortmerna/per_sample/{sample}/out/aligned.log",
    log: "results/3_sortmerna/{sample}_sortmerna.log"
    shell:
        "sortmerna "
        "{params.dbfiles} "
        "--reads {input.r1} "
        "--reads {input.r2} "
        "--threads {threads} "
        "--workdir results/3_sortmerna/per_sample/{wildcards.sample} "
        "--idx-dir resources/sortmerna_db &> {log}; "
        "mv {params.results_log} {output} ; "
        "rm -rf results/3_sortmerna/per_sample/{wildcards.sample}"

rule collect_sortmerna_results:
    input:
        results = get_sortmerna_results()
    output:
        out_figure = report("results/3_sortmerna/rRNA_percentages.png", category="Sortmerna analysis"),
        out_table  = report("results/3_sortmerna/rRNA_percentages.txt", category="Sortmerna analysis"),
    params:
        script = srcdir("../Scripts/gather_sortmerna.py")
    shell:
        "python {params.script} -i {input.results}  "
