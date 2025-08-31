# ORACLE-seq

**Optimal Rarefaction Analysis with Confidence Level Estimation**

🧬 **What it does:** Statistically models feature (isoform/gene) discovery in long-read RNA-seq data to optimize your sequencing strategy

📊 **Input:** Your rarefaction curve data (reads vs features discovered)  
📈 **Output:** Statistical models with confidence intervals + sequencing recommendations

## Quick Example

```r
# 1. Install (one time)
source("run_oracle_seq.R")

# 2. Run analysis
run_oracle_seq(
  data_file = "my_rarefaction_data.rds",
  analysis_name = "my_study", 
  output_dir = "results"
)

# 3. Check results in: results/my_study/
```

**What you get:**
- 📊 Statistical models of your discovery curves
- 📈 Extrapolation predictions with confidence intervals  
- 💰 Cost-benefit analysis for additional sequencing
- 📋 Publication-ready reports and plots

---

## Installation

**Requirements:** R (≥4.0.0)

**Auto-install:** ORACLE-seq automatically installs missing dependencies:
- tidyverse, broom, patchwork, minpack.lm, scales, parallel, ggrepel, tidytext

**Setup:**
```bash
git clone https://github.com/your-username/oracle-seq.git
cd oracle-seq
```

---

## Supported Data Types

✅ **Single-cell RNA-seq** (scRNA-seq)  
✅ **Single-nucleus RNA-seq** (snRNA-seq)  
✅ **Bulk RNA-seq**  
✅ **Direct RNA-seq**  

**Platforms:** Oxford Nanopore (ONT), PacBio

---

## Input Data Format

**File type:** RDS file with a data frame

**Required columns:**
- `curve_id`: Experiment/condition identifier  
- `total_reads_unified`: Total sequencing reads (x-axis)
- `featureNum`: Features discovered at that depth (y-axis)

### Realistic Example Input

```r
# Example: Longbench-style multi-protocol comparison
# Each curve needs multiple points (10-20+ recommended) for good model fitting

rarefaction_data <- data.frame(
  curve_id = rep(c(
    "ONT_Bulk_Genes_Discovery",
    "ONT_SC_Genes_Discovery", 
    "PacBio_Bulk_Genes_Discovery",
    "PacBio_SC_Genes_Discovery"
  ), each = 15),
  
  total_reads_unified = rep(c(
    # Rarefaction depths (15 points per curve)
    50000, 100000, 250000, 500000, 750000, 
    1000000, 1500000, 2000000, 3000000, 4000000,
    5000000, 7500000, 10000000, 15000000, 20000000
  ), 4),
  
  featureNum = c(
    # ONT Bulk Genes (typical saturation curve)
    8500, 12000, 16500, 19000, 20500,
    21500, 22000, 22300, 22600, 22800,
    22900, 22980, 23000, 23020, 23030,
    
    # ONT SC Genes (lower saturation due to dropout)
    5000, 8000, 12000, 15000, 17000,
    18500, 19500, 20000, 20300, 20500,
    20650, 20750, 20800, 20830, 20850,
    
    # PacBio Bulk Genes (higher saturation, cleaner data)
    9000, 13000, 18000, 21000, 23000,
    24500, 25500, 26000, 26400, 26600,
    26750, 26850, 26900, 26930, 26950,
    
    # PacBio SC Genes (better than ONT SC)
    6000, 9500, 14000, 17500, 19500,
    21000, 22000, 22500, 22800, 22950,
    23050, 23120, 23150, 23170, 23180
  )
)

# Save for ORACLE-seq
saveRDS(rarefaction_data, "longbench_rarefaction_data.rds")

# Check structure
str(rarefaction_data)
head(rarefaction_data, 10)
```

### Quick Data Generation Template

```r
# Template for creating your own rarefaction data
create_rarefaction_example <- function() {
  
  # Define your experimental conditions
  conditions <- c("Protocol_A_Genes", "Protocol_B_Genes", "Protocol_A_Isoforms")
  
  # Define rarefaction depths (customize as needed)
  depths <- c(1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7)
  
  # Create empty data frame
  data_list <- list()
  
  for (condition in conditions) {
    for (i in seq_along(depths)) {
      data_list[[length(data_list) + 1]] <- data.frame(
        curve_id = condition,
        total_reads_unified = depths[i],
        featureNum = YOUR_RAREFACTION_FUNCTION(depths[i], condition)
      )
    }
  }
  
  rarefaction_data <- do.call(rbind, data_list)
  saveRDS(rarefaction_data, "my_rarefaction_data.rds")
  return(rarefaction_data)
}
```

**Minimum requirements:**
- ≥3 data points per curve_id (10-20+ recommended)
- Non-negative values
- Generally increasing trend (more reads → more features)

---

## Core Features

### 🎯 **Smart Model Selection**
- **7 competing models:** Linear, Michaelis-Menten, Asymptotic Exponential, Power Law, Logarithmic, Shifted Logarithmic, Hill
- **Automatic best-fit selection** via AIC comparison
- **Bayesian Model Averaging** for robust predictions

### 📊 **Statistical Rigor**  
- **1000+ bootstrap iterations** for confidence intervals
- **Model stability assessment** with confidence categorization:
  - 🟢 **High confidence:** <5% CV  
  - 🟡 **Moderate confidence:** <15% CV
  - 🟠 **Low confidence:** <30% CV
  - 🔴 **Very low confidence:** ≥30% CV

### 💡 **Practical Insights**
- **Marginal returns analysis:** How many new features per million additional reads?
- **Optimal sequencing depth** recommendations
- **Cost-benefit analysis** for experimental design

---

## Example Use Cases

### 🔬 **Data Quality Assessment**
*"Are my features genuine biology or technical artifacts?"*

Create separate rarefaction curves for annotated artifacts vs genuine events, then model their discovery rates to quantify data quality.

### 💰 **Experimental Design**  
*"Is sequencing deeper worth the cost?"*

Model feature discovery to predict gains from additional sequencing and optimize your budget allocation.

### 🧪 **Protocol Comparison**
*"Which sequencing protocol discovers features most efficiently?"*

Compare discovery curves across different protocols, platforms, or experimental conditions.

To Do: ### 📈 **Flexible Discovery Modeling**
Check if model discovery is accurate with feature annotation types like: **
- Differentially expressed genes
- Novel isoforms  
- Intron retention events
- Splice junctions
- Cell type markers

---

## Advanced Usage

### Custom Parameters
```r
run_oracle_seq(
  data_file = "data.rds",
  analysis_name = "detailed_study",
  thresholds = c(20, 10, 5, 1, 0.5, 0.1),    # Custom thresholds
  bootstrap_samples = 2000,                    # More iterations
  confidence_level = 0.99                     # Higher confidence
)
```

### Multiple Analyses
```r
# Wrapper handles name conflicts automatically
run_oracle_seq(data_file = "protocol_A.rds", analysis_name = "comparison")
run_oracle_seq(data_file = "protocol_B.rds", analysis_name = "comparison")  # Creates "comparison_2"
```

---

## Output Structure

```
results/
├── my_study/
│   ├── analysis_config.yaml      # Parameters used
│   ├── statistical_analysis/     # Model results & reports  
│   ├── plots/                   # All visualizations
│   ├── reports/                 # Text summaries
│   └── logs/                    # Analysis logs
```

---

## How It Works

**The ORACLE-seq Approach:**

1. **Input:** Your rarefaction curve coordinates (reads vs features)
2. **Model Competition:** 7 mathematical models compete to best describe your data
3. **Uncertainty Quantification:** Bootstrap resampling generates confidence intervals
4. **Extrapolation:** Predict feature discovery at higher sequencing depths
5. **Recommendations:** Statistical backing for sequencing investment decisions

**Statistical Foundation:**
- Model selection via Akaike Information Criterion (AIC)
- Bayesian Model Averaging accounts for model uncertainty
- Bootstrap confidence intervals (1000+ iterations)
- Coefficient of variation quantifies prediction reliability

---

## Citation

*[Publication details coming soon]*

**Developed for the Longbench project** - benchmarking matched ONT/PacBio datasets across scRNA-seq, snRNA-seq, bulk RNA-seq, and direct RNA-seq from lung carcinoma cell lines.

---

## License

*[License information to be added]*

---

## Contributing

Found a bug? Want a feature? Open an issue or submit a pull request!

**Contact:** mschauhan@student.unimelb.edu.au
