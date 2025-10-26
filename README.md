# Meta-analysis Pipeline at Strain Levels from Shotgun MetaG 
[![Snakemake](https://img.shields.io/badge/workflow-Snakemake-3D9970?style=flat&logo=python&logoColor=white)](https://snakemake.github.io/)
[![FastQC](https://img.shields.io/badge/QC-FastQC-1E90FF?style=flat&logo=data:image/png;base64,iVBORw0KGgo=)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[![MultiQC](https://img.shields.io/badge/QC-MultiQC-FF69B4?style=flat&logo=python&logoColor=white)](https://multiqc.info/)
[![Kraken2](https://img.shields.io/badge/classifier-Kraken2-8A2BE2?style=flat&logo=linux&logoColor=white)](https://ccb.jhu.edu/software/kraken2/)
[![Bracken](https://img.shields.io/badge/classifier-Bracken-FF4500?style=flat&logo=linux&logoColor=white)](https://ccb.jhu.edu/software/bracken/)
[![MIDAS2](https://img.shields.io/badge/metagenomics-MIDAS2-6A5ACD?style=flat&logo=github&logoColor=white)](https://github.com/czbiohub/MIDAS2)

## Project Description:
The project aims to build a reproducible metagenomic pipeline to identify strain-specific changes in type 2 diabetes gut microbiome and metformin treatment.

## Research Questions:
Here we define strain changes as either 
1) Single nucleotide polymorphism (SNPs) within gene of a species AND/OR
2) Changes in gene presence and absence (aka function gain or loss) via copy number variants (CNV) of that gene in a species.
With that goal in mind, we build a workflow with two main steps:
- The first step is to identify differential abundance changes of all species in diabetic condition and treatment.
- With those changed species in mind, the second step is to see if these species also have strains that vary among different sample groups based on the definition mentioned above.

```mermaid
flowchart LR
    A[Step 1 - FastQC] --> B[MultiQC]
    B --> C[BBDuk: Trimming based on Phred Score]
    C --> D[Kraken: Taxa classifiying]
    D --> E[Bracken: Estimating taxa abundance]
    E --> F[phylogenize.v1: Identify increased gene families]
    C --> F[Step 2 - MIDAS2: ID strain variation within each species]
```
## Result Figure:
After species differential abundance analysis, we identify a list of significant species that we cross-compared with literature (upper figure), before they are grouped into phylum level (lower):

<div style="display: flex; justify-content: space-between; align-items: center;">
  <img src="https://github.com/user-attachments/assets/1662cc1f-fe4f-4f68-ac0f-a0b1d1748f79" alt="Screenshot 2025-10-25 at 8 03 13 PM" width="50%">
  <img src="https://github.com/user-attachments/assets/b0f4683e-f21a-4d61-aed4-35ca912c3afa" alt="Screenshot 2025-10-25 at 7 56 50 PM" width="45%">
</div>


We then run a home-built tool -phylogenize- (https://github.com/pbradleylab/phylogenize) to reveal functional changes at species level, corrected for their phylogenetic relationship:
<img width="374" height="566" alt="Screenshot 2025-10-25 at 8 03 37 PM" src="https://github.com/user-attachments/assets/72b98d13-577b-436c-aadf-b1d3506606b6" />

## Conclusion:

## Citations:
