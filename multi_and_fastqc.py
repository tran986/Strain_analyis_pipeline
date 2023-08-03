(DBT,SMPS,PAIR)=glob_wildcards("data/raw/{diabetesType}/{samplename}_{pair}.fastq.gz")

def get_multiqc_subsamples(wildcards):
    metricsLST = []
    for (d,s,p) in zip(DBT,SMPS,PAIR):
      metricsLST.append(rules.run_fastqc.output[0].format(d=d,s=s,p=p))
    return metricsLST

   
#expand("data/processed/fastqc/{d}/{d}_{s}_{p}_fastq_fastqc.html", d=DBT,s=SMPS,p=PAIR),
rule all:
  input: expand("data/processed/multiqc/{d}/multiqc_report.html", d=DBT,s=SMPS,p=PAIR)
     
   
rule rename:
  input:
   "data/raw/{d}/{s}_{p}.fastq.gz"
  output:
   "data/processed/rename/{d}/{d}_{s}_{p}_fastq.gz"
  shell:
   "ln -s $PWD/{input} $PWD/{output}"

rule all_fastqc:
  input:
   [f"data/processed/fastqc/{d}_{s}_{p}_fastqc.html" for (d,s,p) in zip(DBT,SMPS,PAIR)]
    
rule run_fastqc:
   input:
    "data/processed/rename/{d}/{d}_{s}_{p}_fastq.gz"
   output:
    "data/processed/fastqc/{d}/{d}_{s}_{p}_fastq_fastqc.html",
    "data/processed/fastqc/{d}/{d}_{s}_{p}_fastq_fastqc.zip",
   params:
    outdir="data/processed/fastqc/{d}"
   shell:
    "/users/PAS1117/osu9664/eMicro-Apps/FastQC-0.11.8.sif {input} -o {params.outdir}"

 #  [ancient(f"data/processed/fastqc/{samplename}/{d}_{s}_{p}.fastqc.html") for (d,s,p) in zip(DBT,SMPS,PAIR)]
 #[(f"data/processed/fastqc/{d}/{d}_{s}_{p}_fastq_fastqc.html") for (d,s,p) in zip(DBT,SMPS,PAIR)]

rule multiqc:
   input:
    get_multiqc_subsamples
   output:
    "data/processed/multiqc/{d}/multiqc_report.html",
    outdir="data/processed/multiqc/{d}"
   log:
    "logs/{d}.log"
   params:
    outdir="data/processed/multiqc/{d}",
    d="{d}"
   shell:
    "/users/PAS1117/osu9664/eMicro-Apps/MultiQC-1.7.sif data/processed/fastqc/{params.d} -o {output.outdir} 2> {log}"



#command runs on the terminal: $snakemake --profile config multiqc --slurm
#[f"data/processed/fastqc/{d}_{s}_{p}.fastqc.html" for (d,s,p) in zip(DBT,SMP,PAIR)]