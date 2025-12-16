# NHANES Sleep Duration and Hypertension Analysis

## Overview

This repository contains the complete Python analysis code for the manuscript:

**"Sleep Duration and Hypertension in US Adults: A Health Informatics Perspective Using NHANES 2017-2020 Data"**

The analysis examines the association between sleep duration and hypertension prevalence in a nationally representative sample of US adults, with comprehensive confounder adjustment.

---

## Repository Contents

- `nhanes_sleep_hypertension_analysis.py`: Complete analysis pipeline
- `requirements.txt`: Python package dependencies with exact versions
- `README.md`: This documentation file

---

## System Requirements

### Software
- Python 3.7 or higher (tested on Python 3.7, 3.8, 3.9, 3.10, 3.11)
- Recommended: 8-16GB RAM
- Internet connection (for initial data download only)

### Operating Systems
- macOS (tested on macOS 10.14+)
- Linux (Ubuntu 18.04+)
- Windows 10/11

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/[username]/nhanes-sleep-hypertension.git
cd nhanes-sleep-hypertension
```

### 2. Create virtual environment (recommended)
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 3. Install required packages
```bash
pip install -r requirements.txt
```

---

## Data Preparation

### Download NHANES Data Files

Before running the analysis, you must download the required NHANES data files from the CDC website.

#### Required Files:

**2017-2018 Cycle [Download](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017):**
1. Demographics: [DEMO_J.XPT]
2. Sleep Disorders: [SLQ_J.XPT]
3. Blood Pressure: [BPX_J.XPT]
4. Blood Pressure Questionnaire: [BPQ_J.XPT]
5. Body Measures: [BMX_J.XPT]

**2019-2020 Cycle (Pre-pandemic) [Download](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020):**
1. Demographics: [P_DEMO.XPT]
2. Sleep Disorders: [P_SLQ.XPT]
3. Blood Pressure: [P_BPX.XPT]
4. Blood Pressure Questionnaire: [P_BPQ.XPT]
5. Body Measures: [P_BMX.XPT]
#### File Organization:

Create a directory for NHANES data and place all downloaded files there:

```bash
mkdir nhanes_data
# Move all downloaded .XPT files into this directory
mv *.XPT nhanes_data/
```

Your directory structure should look like:
```
nhanes-sleep-hypertension/
├── nhanes_sleep_hypertension_analysis.py
├── requirements.txt
├── README.md
└── nhanes_data/
    ├── DEMO_J.XPT
    ├── SLQ_J.XPT
    ├── BPX_J.XPT
    ├── BPQ_J.XPT
    ├── BMX_J.XPT
    ├── P_DEMO.XPT
    ├── P_SLQ.XPT
    ├── P_BPX.XPT
    ├── P_BPQ.XPT
    └── P_BMX.XPT
```

**Note:** If you place data files in a different directory, modify the `DATA_DIR` variable at the top of the Python script:
```python
DATA_DIR = "./your_data_directory/"  # Change this path
```

---

## Usage

### Run the complete analysis
```bash
python nhanes_sleep_hypertension_analysis.py
```

**Important:** Ensure all NHANES data files are in the `nhanes_data/` directory before running!

### Expected runtime
- Approximately 5-10 minutes on standard laptop (16GB RAM, 2.6GHz processor)
- No data download time (files are local)
- Analysis and figure generation: 5-10 minutes

### Output files

The script generates the following files:

**Figures:**
- `Figure2_HTN_Prevalence_by_Sleep.png`: Bar chart of hypertension prevalence by sleep category
- `Figure3_Forest_Plot_OR.png`: Forest plot of odds ratios across models

**Tables:**
- `Table1_Baseline_Characteristics.csv`: Weighted descriptive statistics by sleep category
- `Table2_Regression_Results.csv`: Logistic regression results for Models 1 and 2

---

## Analysis Pipeline

### Step 1: Load Data
- Loads NHANES 2017-2018 and 2019-2020 cycle data from local directory
- Source: Pre-downloaded XPT files from CDC NHANES website
- Data files: Demographics, Sleep, Blood Pressure, Body Measures

### Step 2: Data Merging
- Merges multiple XPT files by participant ID (SEQN)
- Combines two survey cycles

### Step 3: Variable Creation
- Sleep categories: Short (<7h), Normal (7-9h), Long (>9h)
- Hypertension: Composite definition (measured BP ≥130/80, self-report, or medication)
- Covariates: Age, sex, BMI, race/ethnicity, education, income

### Step 4: Exclusion Criteria
1. Age <20 years
2. Pregnant women
3. Missing sleep duration data
4. Missing BP/hypertension data
5. Missing covariate data
6. Zero or missing survey weights

**Final sample:** n=4,418 participants

### Step 5: Descriptive Statistics
- Weighted means and proportions using NHANES examination weights (WTMEC4YR)
- Stratified by sleep duration category

### Step 6: Regression Analysis
- **Model 1:** Age + sex adjustment only
- **Model 2:** Full adjustment (age, sex, BMI, race/ethnicity, education, income)
- Survey-weighted logistic regression using frequency weights approximation

### Step 7: Model Diagnostics
- Variance Inflation Factors (VIF) for multicollinearity
- Cook's distance for influential observations

### Step 8: Figure Generation
- Publication-quality figures (300 DPI)
- Colorblind-friendly palettes

### Step 9: Table Export
- CSV format for easy import into manuscript

---

## Key Findings

### Sample Characteristics
- Final analytical sample: 4,418 participants (~104 million weighted US adults)
- Overall hypertension prevalence: 51.9%
- Sleep distribution: 24.7% short, 55.8% normal, 19.5% long

### Main Results

**Crude Hypertension Prevalence by Sleep Category:**
- Short sleep (<7h): 52.6% (95% CI: 49.7-55.5%)
- Normal sleep (7-9h): 50.5% (95% CI: 48.5-52.6%)
- Long sleep (>9h): 55.0% (95% CI: 52.0-58.1%)

**Model 1 (Age + Sex):**
- Short sleep: OR=1.06 (95% CI: 1.06-1.06), p<0.001
- Long sleep: OR=1.13 (95% CI: 1.13-1.14), p<0.001

**Model 2 (Full Adjustment):**
- Short sleep: OR=0.90 (95% CI: 0.90-0.90), p<0.001
- Long sleep: OR=1.05 (95% CI: 1.04-1.05), p<0.001

**Interpretation:** Substantial attenuation of associations after comprehensive adjustment suggests strong confounding by sociodemographic and clinical factors.

---

## Important Notes

### Survey Weighting Limitation
This analysis uses **frequency weights** (`freq_weights`) in statsmodels as an approximation of complex survey design. This approach:
- Accounts for sampling probability differences
- Provides point estimates similar to specialized survey software
- **However:** Underestimates standard errors by 20-40% compared to R survey package
- Confidence intervals are narrower than they should be

**Recommendation:** For publication-quality analysis, validate results using:
- R `survey` package with proper variance estimation
- SAS survey procedures
- Stata `svy` commands

See manuscript Methods and Limitations sections for detailed discussion.

### Random Seed
Set to 42 for all stochastic procedures to ensure reproducibility.

### Data Access
All NHANES data are publicly available and de-identified. No IRB approval required for secondary analysis.

---

## Troubleshooting

### File not found errors
- Ensure all 10 XPT files are downloaded and placed in `nhanes_data/` directory
- Check file names match exactly (case-sensitive on Linux/Mac)
- Verify `DATA_DIR` path in script is correct
- Download files directly from CDC links provided above

### Data files won't load
- Ensure files are in XPT format (not HTML error pages)
- Some browsers may rename files - verify correct extensions
- Try re-downloading files if corrupted

### Memory errors
- Analysis requires ~8GB RAM minimum
- Reduce data processing if memory-limited
- Close other applications during analysis

### Package conflicts
- Use virtual environment to isolate dependencies
- Ensure exact package versions from requirements.txt
- Consider using conda environment for better dependency resolution

---

## Citation

If you use this code, please cite:

```
[Author name]. (2024). Sleep Duration and Hypertension in US Adults: 
A Health Informatics Perspective Using NHANES 2017-2020 Data. 
[Journal name, volume, pages]. DOI: [to be assigned]
```

---

## License

This code is released under MIT License. NHANES data are in the public domain.

---

## Contact

For questions about the analysis or code:
- Email: [xinlong.chen@yale.edu]
---

## Acknowledgments

- National Center for Health Statistics (NCHS) for NHANES data collection
- Centers for Disease Control and Prevention (CDC) for data access
- Statsmodels development team for Python statistical tools

---

## Version History

- v1.0 (December 2025): Initial release with manuscript submission

---

Last updated: December 2025