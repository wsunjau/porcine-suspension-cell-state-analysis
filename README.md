# Porcine suspension cell-state analysis
This repository contains the code used to analyse longitudinal single-cell RNA sequencing data from porcine muscle-derived cells during adaptation to serum-free suspension culture.

## Study overview
A satellite-cell-enriched primary culture was progressively adapted from adherent culture to serum-free suspension culture without genetic modification. Single-cell RNA sequencing was performed at nine stages: P5, P7, P9, SD7, SD17, SD38, SD147, SD277 and SD300.
The analyses were used to examine changes in cell composition and transcriptional state during adaptation. The main cell states identified in the study included myogenic, mesenchymal stromal-like, extracellular matrix-remodelling and adipogenic-like populations.

## Repository structure

```text
porcine-suspension-cell-state-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   ├── 01_quality_control.R
│   ├── 02_normalization_integration.R
│   ├── 03_clustering_annotation.R
│   ├── 04_marker_analysis.R
│   ├── 05_differential_expression.R
│   ├── 06_pathway_analysis.R
│   └── 07_pseudotime_analysis.R
└── environment/
    └── sessionInfo.txt
```

## Analysis workflow

The repository includes code for:
1. Quality control and cell filtering
2. Normalization and identification of variable features
3. Data integration across culture stages
4. Dimensionality reduction and clustering
5. Identification and annotation of cell states
6. Marker-gene and differential-expression analyses
7. Functional enrichment analysis
8. Pseudotime analysis



## Reproducibility
To improve reproducibility:
* Relative file paths are used wherever possible.
* Analysis parameters are recorded in the corresponding scripts.
* R and package versions are provided in `environment/sessionInfo.txt`.
* Random seeds are specified for analyses involving stochastic procedures where applicable.


## License
The code in this repository is available under the [MIT License](LICENSE).

## Contact
For questions about the code or analysis, please open an issue in this repository.
