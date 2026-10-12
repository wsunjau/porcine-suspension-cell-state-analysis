# Porcine suspension cell-state analysis
This repository contains the R code used for longitudinal single-cell RNA sequencing analysis of porcine muscle-derived cells during adaptation to serum-free suspension culture.

## Study overview
A satellite-cell-enriched porcine primary culture was progressively adapted from adherent culture to serum-free suspension culture without genetic modification. Single-cell RNA sequencing was performed at nine culture stages: P5, P7, P9, SD7, SD17, SD38, SD147, SD277 and SD300. The analysis examined changes in cell composition and transcriptional state during culture adaptation. The identified cell states included myogenic, mesenchymal stromal-like, extracellular matrix-remodelling and adipogenic-like populations.

## Repository structure

```text
porcine-suspension-cell-state-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── single_cell_analysis.R
└── sessionInfo.txt
```

## Analysis workflow
The complete analysis is provided in a single R script:
```text
single_cell_analysis.R
```

The script is divided into ten sections.

### 01. Packages and settings
Loads the required R packages and defines the settings, colour palettes and parameters used in the analysis.

### 02. Data loading
Imports the single-cell gene-expression data or processed Seurat objects.

### 03. Quality control
Calculates quality-control metrics and filters cells according to the number of detected genes, total transcript counts and mitochondrial transcript proportion.

### 04. Normalization and scaling
Performs data normalization, variable-feature selection, data integration, scaling and principal component analysis.

### 05. Dimensionality reduction and clustering
Constructs the nearest-neighbour graph, identifies cell clusters and generates UMAP representations.

### 06. Marker analysis
Identifies marker genes for individual clusters and generates marker-expression plots.

### 07. Cell-state annotation
Assigns cell-state identities to clusters based on marker-gene expression and transcriptional characteristics.

### 08. Pathway analysis
Performs functional enrichment analysis using genes identified from the relevant cell-state or differential-expression comparisons.

### 09. Cell-cycle analysis
Scores cell-cycle-associated gene expression and compares cell-cycle distributions across culture stages and cell states.

### 10. Pseudotime analysis
Reconstructs cell-state trajectories and examines transcriptional changes along pseudotime.


## Running the analysis
Clone or download the repository:
```bash
git clone https://github.com/Leonrunning/porcine-suspension-cell-state-analysis.git
cd porcine-suspension-cell-state-analysis
```

Open the following file in RStudio:

```text
single_cell_analysis.R
```
The input and output paths in the script should be changed to match the corresponding local file locations.


## Software environment

The analyses were performed in R. Information on the R version, operating system and package versions used for the final analysis is provided in:
```text
sessionInfo.txt
```

Differences in software or package versions may result in minor differences in the output.

## Reproducibility
The script records the principal parameters used for data processing and downstream analysis. Relative file paths are used where possible, and random seeds are specified for procedures involving stochastic calculations.

Before running the script, users should confirm that:
- The required R packages are installed.
- Input file names and paths are correct.
- The required output directories have been created.
- Sufficient memory is available for processing the single-cell objects.

## Citation
If you use the code in this repository, please cite the associated paper:
```text
Citation information will be added after publication.
```

## License
The code in this repository is available under the [MIT License](LICENSE).

## Contact
For questions about the code or analysis, please open an issue in this repository.
