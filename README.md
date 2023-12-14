# RayFin-TEvolution

# Workflow Documentation for "Swimming to Jumping: Unveiling the Dynamic World of Transposable Elements in Ray-finned Fish Evolution"

## Overview
These are the steps followed to obtain the results presented in the paper.

### Initial Data Processing
- The data was provided in the format .fsa.
- Due to contaminations like adapters and duplicate sequences, the file was initially rejected by NCBI.

### Genome Annotation
- **Script Used**: `fullbrakerpipeline_bichir.sh`
- **Purpose**: To annotate the genome.
- **Outcome**: Part of this process involves running RepeatModeler, which generates two key files: `.tbl` and `.out`. The `.tbl` file contains the percentage of different repeats in the genome, while the `.out` file details the position, length, and other aspects of the repeats.

### Transposable Elements (TEs) Analysis
- **Script**: `getTEdetails.py`
  - **Input**: `.tbl` file.
  - **Function**: Prints TE details in tabular format.
- **Script**: `TEcontentVSgenomesize.py`
  - **Input**: `.tbl` file and genome size.
  - **Function**: Calculates the percentage of repeats in the genome.
- **Script**: `getTEdetails_everyfamily.py`
  - **Input**: `.out` file.
  - **Function**: Prints the details of each family in tabular format.

### Pfam Analysis and Result Plotting
- **Script Used**: `Pfam-analysis.sh`
- **Purpose**: Conducts Pfam analysis.
- **R Scripts for Plotting**:
  - `boxplot.R`  
  - `LinearModelsForTEs.R`
  - `phylopca_bichir_AD.R` — Creates the phylomorphospace.
  - `TE_types_allfishesplot.R` — Calculates phylogenetic generalized least squares (PGLS) and generates plots.
  - `t-testandlambatest_bichir.R` — Calculates Pagel's lambda and Blomberg's K.
  - `LinearModels_for_TEs_OU.R` — Calculates Ornstein-Uhlenbeck (OU) models.
  - `model_fitting_TEclasses.R` — Implements phylostep model for stepwise PGLS.
  - `TreePrAbHeatmap.R` — Uses SIMMAP code for ancestral state reconstruction.

## Getting Started
Clone this repository and run the scripts in an R or Python environment as applicable. Ensure all necessary dependencies are installed.

## Contributions
Contributions are welcome. Please adhere to coding standards and submit pull requests for review.

---

Contact the repository maintainers for further information or support.
