[![Snakemake](https://img.shields.io/badge/workflow-Snakemake-3D9970?style=flat&logo=python&logoColor=white)](https://snakemake.github.io/)
[![FastQC](https://img.shields.io/badge/QC-FastQC-1E90FF?style=flat&logo=data:image/png;base64,iVBORw0KGgo=)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[![MultiQC](https://img.shields.io/badge/QC-MultiQC-FF69B4?style=flat&logo=python&logoColor=white)](https://multiqc.info/)
[![Kraken2](https://img.shields.io/badge/classifier-Kraken2-8A2BE2?style=flat&logo=linux&logoColor=white)](https://ccb.jhu.edu/software/kraken2/)
[![Bracken](https://img.shields.io/badge/classifier-Bracken-FF4500?style=flat&logo=linux&logoColor=white)](https://ccb.jhu.edu/software/bracken/)
[![MIDAS2](https://img.shields.io/badge/metagenomics-MIDAS2-6A5ACD?style=flat&logo=github&logoColor=white)](https://github.com/czbiohub/MIDAS2)

```mermaid
flowchart LR
    A[FastQC] --> B[MultiQC]
    B --> C[BBDuk: Trimming based on Phred Score]
    C --> D[Kraken: Taxa classifiying]
    D --> E[Bracken: Estimating taxa abundance]
    C --> F[MIDAS2: ID strain variation within each species]
