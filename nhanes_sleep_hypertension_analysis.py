"""
NHANES Sleep Duration and Hypertension Analysis
Complete analysis pipeline for manuscript submission

Author: Xinlong Chen
Date: December 2025
Python Version: 3.7+

This script performs the complete analysis reported in the manuscript:
"Sleep Duration and Hypertension in US Adults: A Health Informatics 
Perspective Using NHANES 2017-2020 Data"

Analysis Steps:
1. Data loading from local NHANES files
2. Data merging and cleaning
3. Exclusion criteria application
4. Weighted descriptive statistics
5. Survey-weighted logistic regression
6. Model diagnostics
7. Figure generation
8. Table output

Note: Survey weighting uses frequency weights approximation.
See manuscript Methods section for detailed discussion of limitations.
"""
# =============================================================================
# 1) Imports & setup
# =============================================================================

import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings("ignore")
np.random.seed(42)

pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)

print("=" * 80)
print("NHANES Sleep Duration and Hypertension Analysis")
print("=" * 80)
print(f"Python: {sys.version.split()[0]}")
print(f"Pandas: {pd.__version__}")
print(f"NumPy: {np.__version__}")
print(f"Statsmodels: {sm.__version__}")

# =============================================================================
# 2) Load NHANES data files (.XPT)
#    Update DATA_PATH to your local folder containing downloaded files
# =============================================================================

DATA_PATH = "./nhanes_data/"  # e.g., "/Users/you/nhanes_data/"

# 2017–2018
demo_2017 = pd.read_sas(f"{DATA_PATH}DEMO_J.xpt")
slq_2017  = pd.read_sas(f"{DATA_PATH}SLQ_J.xpt")
bpq_2017  = pd.read_sas(f"{DATA_PATH}BPQ_J.xpt")
bpx_2017  = pd.read_sas(f"{DATA_PATH}BPX_J.xpt")
bmx_2017  = pd.read_sas(f"{DATA_PATH}BMX_J.xpt")

# 2019–2020 pre-pandemic
demo_2019 = pd.read_sas(f"{DATA_PATH}P_DEMO.xpt")
slq_2019  = pd.read_sas(f"{DATA_PATH}P_SLQ.xpt")
bpq_2019  = pd.read_sas(f"{DATA_PATH}P_BPQ.xpt")
bpx_2019  = pd.read_sas(f"{DATA_PATH}P_BPXO.xpt")  # oscillometric BP file
bmx_2019  = pd.read_sas(f"{DATA_PATH}P_BMX.xpt")

print("\nLoaded files successfully.")

# =============================================================================
# 3) Merge within each cycle, then combine cycles
# =============================================================================

df_2017 = (
    demo_2017
    .merge(slq_2017, on="SEQN", how="left")
    .merge(bpq_2017, on="SEQN", how="left")
    .merge(bpx_2017, on="SEQN", how="left")
    .merge(bmx_2017, on="SEQN", how="left")
)
df_2017["CYCLE"] = "2017-2018"

df_2019 = (
    demo_2019
    .merge(slq_2019, on="SEQN", how="left")
    .merge(bpq_2019, on="SEQN", how="left")
    .merge(bpx_2019, on="SEQN", how="left")
    .merge(bmx_2019, on="SEQN", how="left")
)
df_2019["CYCLE"] = "2019-2020"

df = pd.concat([df_2017, df_2019], ignore_index=True)
print(f"\nCombined dataset: {df.shape[0]:,} rows, {df.shape[1]:,} columns")

# =============================================================================
# 4) Adjust weights for combined 4-year period
#    NHANES guidance: WTMEC4YR = WTMEC2YR * (2/4)
# =============================================================================

df["WTMEC4YR"] = df["WTMEC2YR"] * (2 / 4)

# =============================================================================
# 5) Create analysis variables
# =============================================================================

# Sleep hours (NHANES naming differs by cycle; keep a unified field)
if "SLD012" in df.columns:
    df["sleep_hours"] = df["SLD012"]
elif "SLQ300" in df.columns:
    df["sleep_hours"] = df["SLQ300"]
else:
    raise KeyError("Sleep duration variable not found (expected SLD012 or SLQ300).")

df["sleep_category"] = pd.cut(
    df["sleep_hours"],
    bins=[0, 7, 9, 24],
    labels=["Short (<7h)", "Normal (7-9h)", "Long (>9h)"],
    right=False
)

# Hypertension composite:
# (1) self-report BPQ020==1 OR (2) measured mean BP >=130/80 OR (3) BP meds BPQ050A==1
df["htn_self_report"] = (df["BPQ020"] == 1).astype(float)
df["htn_medication"]  = (df["BPQ050A"] == 1).astype(float)

# Measured BP: standard (2017–2018) + oscillometric (2019–2020) fallback
df["avg_systolic"] = np.nan
df["avg_diastolic"] = np.nan

if {"BPXSY1", "BPXSY2", "BPXSY3"}.issubset(df.columns):
    df["avg_systolic"] = df[["BPXSY1", "BPXSY2", "BPXSY3"]].mean(axis=1)
if {"BPXDI1", "BPXDI2", "BPXDI3"}.issubset(df.columns):
    df["avg_diastolic"] = df[["BPXDI1", "BPXDI2", "BPXDI3"]].mean(axis=1)

if {"BPXOSY1", "BPXOSY2", "BPXOSY3"}.issubset(df.columns):
    mask = df["avg_systolic"].isna()
    df.loc[mask, "avg_systolic"] = df.loc[mask, ["BPXOSY1", "BPXOSY2", "BPXOSY3"]].mean(axis=1)

if {"BPXODI1", "BPXODI2", "BPXODI3"}.issubset(df.columns):
    mask = df["avg_diastolic"].isna()
    df.loc[mask, "avg_diastolic"] = df.loc[mask, ["BPXODI1", "BPXODI2", "BPXODI3"]].mean(axis=1)

df["htn_measured"] = ((df["avg_systolic"] >= 130) | (df["avg_diastolic"] >= 80)).astype(float)
df["hypertension"] = (
    (df["htn_self_report"] == 1) |
    (df["htn_measured"] == 1) |
    (df["htn_medication"] == 1)
).astype(float)

# Core covariates
df["age"] = df["RIDAGEYR"]
df["female"] = (df["RIAGENDR"] == 2).astype(int)
df["bmi"] = df["BMXBMI"]

df["race_ethnicity"] = df["RIDRETH3"].map({
    1: "Hispanic",
    2: "Hispanic",
    3: "Non-Hispanic White",
    4: "Non-Hispanic Black",
    6: "Non-Hispanic Asian",
    7: "Other/Multiracial",
})

df["education"] = df["DMDEDUC2"].map({
    1: "Less than HS",
    2: "Less than HS",
    3: "HS/GED",
    4: "Some college",
    5: "College grad+",
})

df["pir"] = df["INDFMPIR"]

# =============================================================================
# 6) Apply exclusions
# =============================================================================

n_initial = len(df)

# Adults only (you used >=20; keep consistent with DMDEDUC2 availability)
df = df[df["age"] >= 20].copy()
n_after_age = len(df)

# Pregnancy exclusion (if available)
if "RIDEXPRG" in df.columns:
    df = df[(df["RIDEXPRG"] != 1) | (df["RIDEXPRG"].isna())].copy()
n_after_preg = len(df)

# Drop missing outcome/exposure/covariates
df = df[df["sleep_hours"].notna()].copy()
n_after_sleep = len(df)

df = df[df["hypertension"].notna()].copy()
n_after_bp = len(df)

df = df[
    df["age"].notna() &
    df["female"].notna() &
    df["bmi"].notna() &
    df["race_ethnicity"].notna() &
    df["education"].notna() &
    df["pir"].notna()
].copy()
n_final = len(df)

# Weight filtering (critical for weighted analysis)
df = df[(df["WTMEC4YR"] > 0) & (df["WTMEC4YR"].notna())].copy()

print("\nSample flow:")
print(f"Initial: {n_initial:,}")
print(f"After age>=20: {n_after_age:,}")
print(f"After pregnancy excl: {n_after_preg:,}")
print(f"After sleep non-missing: {n_after_sleep:,}")
print(f"After outcome non-missing: {n_after_bp:,}")
print(f"After covariates complete: {n_final:,}")
print(f"Final after weight filter: {len(df):,}")

# =============================================================================
# 7) Weighted descriptive statistics (approximate SEs)
# =============================================================================

def weighted_mean_se(x, w):
    m = np.average(x, weights=w)
    v = np.average((x - m) ** 2, weights=w)
    se = np.sqrt(v / len(x))
    return m, se

def weighted_proportion(x01, w):
    p = np.average(x01, weights=w)
    se = np.sqrt(p * (1 - p) / len(x01))
    return p, se

age_mean, age_se = weighted_mean_se(df["age"], df["WTMEC4YR"])
htn_p, htn_se = weighted_proportion(df["hypertension"], df["WTMEC4YR"])

print("\nWeighted summary (approx):")
print(f"Age mean (SE): {age_mean:.1f} ({age_se:.2f})")
print(f"Hypertension % (SE): {htn_p*100:.1f}% ({htn_se*100:.2f}%)")

# Prevalence by sleep category
results_by_sleep = []
for cat in ["Short (<7h)", "Normal (7-9h)", "Long (>9h)"]:
    sub = df[df["sleep_category"] == cat]
    if len(sub) == 0:
        continue
    p, se = weighted_proportion(sub["hypertension"], sub["WTMEC4YR"])
    results_by_sleep.append({
        "Category": cat,
        "Prevalence": p * 100,
        "SE": se * 100,
        "CI_Lower": (p - 1.96 * se) * 100,
        "CI_Upper": (p + 1.96 * se) * 100,
        "N": len(sub)
    })
results_df = pd.DataFrame(results_by_sleep)

# =============================================================================
# 8) Survey-weighted logistic regression (freq_weights approximation)
# =============================================================================

df["sleep_short"] = (df["sleep_category"] == "Short (<7h)").astype(int)
df["sleep_long"] = (df["sleep_category"] == "Long (>9h)").astype(int)

df["race_black"] = (df["race_ethnicity"] == "Non-Hispanic Black").astype(int)
df["race_hispanic"] = (df["race_ethnicity"] == "Hispanic").astype(int)
df["race_asian"] = (df["race_ethnicity"] == "Non-Hispanic Asian").astype(int)
df["race_other"] = (df["race_ethnicity"] == "Other/Multiracial").astype(int)

df["edu_less_hs"] = (df["education"] == "Less than HS").astype(int)
df["edu_hs"] = (df["education"] == "HS/GED").astype(int)
df["edu_some_college"] = (df["education"] == "Some college").astype(int)

reg_vars = [
    "hypertension", "sleep_short", "sleep_long", "age", "female", "bmi",
    "race_black", "race_hispanic", "race_asian", "race_other",
    "edu_less_hs", "edu_hs", "edu_some_college", "pir", "WTMEC4YR"
]
df_reg = df[reg_vars].dropna().copy()

y = df_reg["hypertension"]
w = df_reg["WTMEC4YR"]

# Model 1: age + sex
X1 = sm.add_constant(df_reg[["sleep_short", "sleep_long", "age", "female"]])
m1 = sm.GLM(y, X1, family=sm.families.Binomial(), freq_weights=w).fit()

# Model 2: fully adjusted (primary)
X2 = sm.add_constant(df_reg[
    ["sleep_short", "sleep_long", "age", "female", "bmi",
     "race_black", "race_hispanic", "race_asian", "race_other",
     "edu_less_hs", "edu_hs", "edu_some_college", "pir"]
])
m2 = sm.GLM(y, X2, family=sm.families.Binomial(), freq_weights=w).fit()

def or_ci(model, term):
    OR = np.exp(model.params[term])
    ci = np.exp(model.conf_int().loc[term])
    p = model.pvalues[term]
    return OR, ci[0], ci[1], p

or1_s, lo1_s, hi1_s, p1_s = or_ci(m1, "sleep_short")
or1_l, lo1_l, hi1_l, p1_l = or_ci(m1, "sleep_long")
or2_s, lo2_s, hi2_s, p2_s = or_ci(m2, "sleep_short")
or2_l, lo2_l, hi2_l, p2_l = or_ci(m2, "sleep_long")

print("\nWeighted logistic regression (freq_weights approximation):")
print(f"Model 1 short sleep OR: {or1_s:.2f} ({lo1_s:.2f}-{hi1_s:.2f}), p={p1_s:.4f}")
print(f"Model 1 long  sleep OR: {or1_l:.2f} ({lo1_l:.2f}-{hi1_l:.2f}), p={p1_l:.4f}")
print(f"Model 2 short sleep aOR: {or2_s:.2f} ({lo2_s:.2f}-{hi2_s:.2f}), p={p2_s:.4f}")
print(f"Model 2 long  sleep aOR: {or2_l:.2f} ({lo2_l:.2f}-{hi2_l:.2f}), p={p2_l:.4f}")

# =============================================================================
# 9) Diagnostics: VIF for multicollinearity
# =============================================================================

vif = pd.DataFrame({
    "Variable": X2.columns[1:],
    "VIF": [variance_inflation_factor(X2.values, i) for i in range(1, X2.shape[1])]
}).sort_values("VIF", ascending=False)

print("\nVIF (Model 2 predictors):")
print(vif.to_string(index=False))

# =============================================================================
# 10) Figures: Prevalence plot + Forest plot (publication-ready output)
# =============================================================================

sns.set_style("whitegrid")

# Figure: prevalence by sleep category
if len(results_df) > 0:
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(results_df["Category"]))
    yv = results_df["Prevalence"].values
    err = (results_df["SE"].values) * 1.96

    ax.bar(x, yv, yerr=err, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(results_df["Category"])
    ax.set_xlabel("Sleep duration category")
    ax.set_ylabel("Hypertension prevalence (%)")
    ax.set_title("Weighted hypertension prevalence by sleep duration (NHANES 2017–Mar 2020)")
    plt.tight_layout()
    plt.savefig("Figure2_HTN_Prevalence_by_Sleep.png", dpi=300, bbox_inches="tight")

# Forest plot: ORs from Model 1 & Model 2
fig, ax = plt.subplots(figsize=(8, 5))
models = ["Model 1", "Model 2"]
ypos = np.arange(len(models))

short_ors = [or1_s, or2_s]
short_lo = [lo1_s, lo2_s]
short_hi = [hi1_s, hi2_s]

long_ors = [or1_l, or2_l]
long_lo = [lo1_l, lo2_l]
long_hi = [hi1_l, hi2_l]

ax.errorbar(short_ors, ypos - 0.1,
            xerr=[np.array(short_ors) - np.array(short_lo),
                  np.array(short_hi) - np.array(short_ors)],
            fmt="o", capsize=4, label="Short (<7h)")

ax.errorbar(long_ors, ypos + 0.1,
            xerr=[np.array(long_ors) - np.array(long_lo),
                  np.array(long_hi) - np.array(long_ors)],
            fmt="s", capsize=4, label="Long (>9h)")

ax.axvline(1.0, linestyle="--")
ax.set_yticks(ypos)
ax.set_yticklabels(models)
ax.set_xlabel("Odds ratio (95% CI)")
ax.set_title("Sleep duration and hypertension (reference: 7–9h)")
ax.legend()
plt.tight_layout()
plt.savefig("Figure3_Forest_Plot_OR.png", dpi=300, bbox_inches="tight")

# =============================================================================
# 11) Tables & outputs
# =============================================================================

# Table 1: weighted baseline characteristics by sleep category (approx SEs)
table1_rows = []
for cat in ["Short (<7h)", "Normal (7-9h)", "Long (>9h)"]:
    sub = df[df["sleep_category"] == cat]
    if len(sub) == 0:
        continue
    wm = sub["WTMEC4YR"]
    age_m, age_se = weighted_mean_se(sub["age"], wm)
    bmi_m, bmi_se = weighted_mean_se(sub["bmi"], wm)
    female_p, female_se = weighted_proportion(sub["female"], wm)
    htn_p, htn_se = weighted_proportion(sub["hypertension"], wm)
    table1_rows.append({
        "Sleep Category": cat,
        "N": len(sub),
        "Age (mean±SE)": f"{age_m:.1f}±{age_se:.2f}",
        "Female %": f"{female_p*100:.1f}",
        "BMI (mean±SE)": f"{bmi_m:.1f}±{bmi_se:.2f}",
        "HTN %": f"{htn_p*100:.1f}",
    })
table1 = pd.DataFrame(table1_rows)

table2 = pd.DataFrame({
    "Model": ["Model 1", "Model 2 (Primary)"],
    "Short Sleep OR": [f"{or1_s:.2f}", f"{or2_s:.2f}"],
    "Short 95% CI": [f"{lo1_s:.2f}-{hi1_s:.2f}", f"{lo2_s:.2f}-{hi2_s:.2f}"],
    "Short p": [f"{p1_s:.4f}", f"{p2_s:.4f}"],
    "Long Sleep OR": [f"{or1_l:.2f}", f"{or2_l:.2f}"],
    "Long 95% CI": [f"{lo1_l:.2f}-{hi1_l:.2f}", f"{lo2_l:.2f}-{hi2_l:.2f}"],
    "Long p": [f"{p1_l:.4f}", f"{p2_l:.4f}"],
})

table1.to_csv("Table1_Baseline_Characteristics.csv", index=False)
table2.to_csv("Table2_Regression_Results.csv", index=False)
vif.to_csv("VIF_Multicollinearity.csv", index=False)
results_df.to_csv("HTN_Prevalence_by_Sleep.csv", index=False)

with open("Model2_Full_Summary.txt", "w") as f:
    f.write(m2.summary().as_text())

print("\nSaved outputs:")
print(" - Figure2_HTN_Prevalence_by_Sleep.png")
print(" - Figure3_Forest_Plot_OR.png")
print(" - Table1_Baseline_Characteristics.csv")
print(" - Table2_Regression_Results.csv")
print(" - VIF_Multicollinearity.csv")
print(" - HTN_Prevalence_by_Sleep.csv")
print(" - Model2_Full_Summary.txt")

print("\nDone.")
