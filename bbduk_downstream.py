#making a file for snakemake files
#using BBduk for quality control for fastq files
#Type2Diabetes_Forslund
#list of [Dir1, Dir2]
dir = ["1","2","3","MHD"]
dirsamples =[]
for d in dir:
    (SMPS, PAIR)= glob_wildcards("data/raw/diabetesType2_" + d + "/{sample}_{pairs}.fastq.gz")
    print(SMPS)
    print(PAIR)
    dirsamples.extend(["diabetesType2_"+ d + f"/{sample}" for sample in SMPS])

print(list(set(dirsamples)))

rule all_bbduk:
  input:
    expand("data/processed/bbduk_downstream/{DIRS}_{pair}.fastq",
       DIRS=dirsamples,
       pair=[1,2]),

rule bbduk:
 input:
  in1 = "data/raw/{dirsamples}_1.fastq.gz",
  in2 = "data/raw/{dirsamples}_2.fastq.gz",
 output:
  out1 = "data/processed/bbduk_downstream/{dirsamples}_1.fastq",
  out2 = "data/processed/bbduk_downstream/{dirsamples}_2.fastq",
 resources:
  time = "0:5:00",
  nodes = 1,
 shell:
  "/users/PAS1117/osu9664/eMicro-Apps/BBTools-38.69.sif bbduk.sh in={input.in1} in2={input.in2} out={output.out1} out2={output.out2} qtrim=r trimq=10"
rule list_dir:
  output: "data/test/list_dir.txt"
  shell: "ls > {output}"


#code needs to be run on terminal: $ snakemake --profile config all_bbduk --slurm
#running with 1 pair of sample: adding [1] at the end of rule all