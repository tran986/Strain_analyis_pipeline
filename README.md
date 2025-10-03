```mermaid
flowchart LR
    A[FastQC] --> B[MultiQC]
    B --> C[BBDuk]
    C --> D[Kraken]
    D --> E[Bracken]
    C --> F[MIDAS2: Strain Analysis]
