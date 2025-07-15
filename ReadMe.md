# ORACLE-seq: Optimal Rarefaction Analysis with Confidence Level Estimation

## Overview

ORACLE-seq is a statistically rigorous R pipeline for **flexible discovery modelling** and **experimental design optimisation** in long-read RNA-sequencing datasets. By modelling the discovery curves of any filtered feature set from your count matrix, ORACLE-seq helps you make informed decisions about current data quality and future sequencing investments.

## Core Concept

**Input**: x,y coordinates (reads vs. features discovered) from any rarefaction analysis  
**Output**: Statistical models with confidence intervals and optimal sequencing recommendations  
**Power**: What you choose to "count" as a "feature" is completely up to your biological question

The tool's flexibility depends on how creative you are as a bioinformatician to ask different biological questions of interest for you to make decisions on with a statistically rigorous method.

## Key Features

- **7-Model Competition**: Linear, Michaelis-Menten, Asymptotic Exponential, Power Law, Logarithmic, Shifted Logarithmic, and Hill models compete via Bayesian Model Averaging
- **Bootstrap Uncertainty Quantification**: 1000+ iterations with confidence intervals
- **Model Stability Assessment**: Coefficient of variation (CV) from bootstrap predictions categorises extrapolation confidence (High: <5% CV, Moderate: <15% CV, Low: <30% CV, Very Low: ≥30% CV)
- **Configurable Thresholds**: Marginal returns analysis at multiple discovery rates
- **Publication-Ready Outputs**: Comprehensive statistical analysis and visualisations
- **Long-Read Optimised**: Designed for ONT/PacBio scRNA-seq, snRNA-seq, bulk RNA-seq, and direct RNA-seq

## Use Cases & Applications

### Data Quality Assessment
*"How much of your count matrix discovery is dominated by internal priming artifacts vs genuine intron retention events (if that's of interest to your experiment). Create separate rarefaction curves for each annotation category and model their discovery rates to quantify the relative contribution of artifacts vs genuine events"*

ORACLE-seq provides quantitative data quality insights by:
- **Comparative discovery modelling**: Separate statistical models for artifact vs genuine feature discovery curves reveal which dominates at different sequencing depths
- **Plateau analysis**: Determine if artifact discovery saturates early while genuine events continue accumulating, indicating data quality vs biological signal
- **Proportional assessment**: Calculate the statistical ratio of artifact to genuine feature discovery with confidence intervals
- **Depth-dependent quality metrics**: Understand how data quality changes with sequencing depth to optimise filtering strategies

### Experimental Design Optimisation  
*"If I were to sequence a couple million more reads or reach some discovery threshold would that be worth it?"*

ORACLE-seq provides statistically-backed answers by:
- **Extrapolating beyond current data**: Models predict feature discovery at higher sequencing depths with confidence intervals
- **Quantifying marginal returns**: Calculates how many new features you'd discover per additional million reads at configurable thresholds (e.g., 10, 1, 0.1 features/million reads)
- **Assessing prediction reliability**: Model stability indicators tell you how much confidence to place in extrapolation estimates
- **Cost-benefit analysis**: Compare predicted discovery gains against sequencing costs to optimise experimental budgets

### Flexible Discovery Modelling
Model the discovery of any annotated feature set:
- Differentially expressed genes
- Novel isoforms 
- Intron retention events
- Technical artifacts
- Splice junctions
- Cell type markers
- ...or any custom annotation strategy

## Technical Approach

**Statistical Rigour**: Multiple mathematical models compete to best describe your discovery curve, with model selection via AIC and uncertainty quantification through bootstrapping.

**Model Competition**: Seven mathematical models simultaneously fit your data, with Bayesian Model Averaging providing robust predictions that account for model uncertainty.

**Confidence Estimation**: Bootstrap resampling (1000+ iterations) generates confidence intervals around predictions, giving you statistical backing for experimental design decisions.

**Model Stability Proxy**: Coefficient of variation (CV = standard deviation / mean) from bootstrap predictions quantifies how much "stake" you should put in the model's extrapolation at each threshold. Lower CV indicates more reliable extrapolation confidence.

## Target Data Types

- **Single-cell RNA-seq** (scRNA-seq)
- **Single-nucleus RNA-seq** (snRNA-seq) 
- **Bulk RNA-seq**
- **Direct RNA-seq**
- **Long-read sequencing platforms**: Oxford Nanopore (ONT), PacBio

## Integration & Workflow

ORACLE-seq integrates into custom analysis pipelines as a discovery modelling engine. You're limited only by your creativity in annotating transcripts and genes in a count matrix - then modelling the discovery of whatever you've annotated and filtered prior to inputting those x,y rarefied values into ORACLE-seq.

## Quick Start

```r
source("run_oracle_seq.R")

run_oracle_seq(
  data_file = "rarefaction_results.rds",
  analysis_name = "intron_retention_discovery", 
  output_dir = "oracle_results",
  thresholds = c(10, 1, 0.1),
  bootstrap_samples = 1000,
  use_bayesian_model_averaging = TRUE,
  confidence_level = 0.95
)
```

## Input Data Requirements

**File Format**: RDS file containing a data frame

**Required Columns**:
- `curve_id`: Character/factor identifying each experimental condition or sample
- `total_reads_unified`: Numeric vector of total sequencing reads (x-axis values)  
- `featureNum`: Numeric vector of features discovered at each depth (y-axis values)

**Example Data Structure**:
```r
# Example rarefaction data frame
data.frame(
  curve_id = c("bulk_RNA_seq", "bulk_RNA_seq", "scRNA_seq", "scRNA_seq"),
  total_reads_unified = c(1000000, 5000000, 500000, 2000000),
  featureNum = c(15000, 18000, 8000, 12000)
)

# Save as RDS for ORACLE-seq input
saveRDS(rarefaction_data, "my_rarefaction_data.rds")
```

**Data Quality Requirements**:
- Minimum 3 data points per `curve_id` for model fitting
- Non-negative values for both depth and feature counts
- Sufficient range in x-values (depth) for meaningful model fitting
- Generally increasing trend (more reads → more features discovered)

**Preprocessing Notes**:
- Multiple curves can be analyzed simultaneously by using different `curve_id` values
- Feature filtering strategy determines what biological question you're modelling

## Citation

*[Publication details to be added]*

## License

*[License information to be added]* 