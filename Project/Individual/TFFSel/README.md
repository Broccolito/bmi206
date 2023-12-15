# High-Altitude Natural Selection Near Transcription Factor Footprints

## Overview
This repository contains the R code used in the study "High-altitude Natural Selection Near Transcription Factor Footprints" by Wanjun Gu. The study explores genetic adaptations in high-altitude environments, particularly focusing on transcription factor (TF) footprints and their relation to natural selection in Andean highlanders.

## Abstract
The research builds upon Vierstra et al.'s creation of a comprehensive map of TF footprints across the human genome, examining if TF-occupied areas are hotspots for genetic adaptations in high-altitude environments. The study utilizes a selection scan on whole genomes of native Andean highlanders and conducts an enrichment test on these genomes to analyze genetic markers in relation to TF footprints and DNAse I hypersensitivity sites.

## Repository Contents

### R Scripts
- `annotate_overlap.R`: Annotates the overlap between genetic markers and TF footprint regions.
- `cleanup.R`: Cleans and prepares the genomic data for analysis.
- `count_overlap.R`: Counts the overlapping genetic markers in designated genomic regions.
- `enrichment.R`: Performs the enrichment test on the Andean genome markers.
- `statistical_analysis.R`: Conducts statistical analysis of the data, including the calculation of p-values and Bayes factors.
- `visualize.R`: Generates visualizations for the study, such as graphs and network diagrams.

### Figures
- Figure XA: Enrichment in DNase I hypersensitivity sites and TF footprints.
- Figure XB: Biological pathway analysis highlighting DNA regulation, replication, and repair.
- Figure XC: Network diagram of pathways related to DNA regulation, replication, and repair.

## Data Source
The data for this study was primarily derived from the supplementary data provided by Vierstra et al. and native Andean highlanders' genomes available upon request by the Simonson lab at UC San Diego. 

## Usage
Each script is designed for a specific part of the data analysis workflow:
1. **Data Preparation**: Run `cleanup.R` to prepare the genomic data.
2. **Overlap Annotation**: Use `annotate_overlap.R` to identify overlaps between genetic markers and TF regions.
3. **Overlap Counting**: Execute `count_overlap.R` to quantify these overlaps.
4. **Enrichment Analysis**: `enrichment.R` is used for conducting the enrichment test.
5. **Statistical Analysis**: `statistical_analysis.R` provides the statistical backbone for the study.
6. **Visualization**: Lastly, `visualize.R` is used to create visual representations of the data and findings.

## Dependencies
- Cytoscape
- BiNGO
- Other R packages as required by the scripts

## Contact
Wanjun Gu: wanjun.gu@ucsf.edu

