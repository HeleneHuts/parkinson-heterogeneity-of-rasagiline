# parkinson-heterogeneity-of-rasagiline

Code for *Causal machine learning reveals heterogeneity in the effect of rasagiline in Parkinson’s disease* (Huts et al.).

This analysis uses causal forests to estimate the individualized treatment effect of rasagiline on Movement Disorder Society-Unified Parkinson's Disease Rating Scale (MDS-UPDRS) part 2+3 score. The analysis focuses on the change in score from baseline to week 14. Individualized treatment effect is expressed as difference in outcome for a patient if they had received rasagiline compared to if they had received placebo. The analysis searches for heterogeneity in the effect of rasagiline and for baseline characteristics that are associated with an enhanced response to treatment.

## Data

This analysis uses data from the Integrated Parkinson’s Database of the Critical Path Institute. The analysis focuses on data from three randomized controlled trials:
- NCT00256204 ADAGIO, a trial in patients with early Parkinson's disease
- NCT02337725 a Japanese trial in patients with early Parkinson's disease
- NCT02337738 a Japanese trial in patients with advanced Parkinson's disease

**The trial data is not included in this repository and cannot be redistributed.** A data access request can be submitted to the Critical Path Institute via https://c-path.org/tools-platforms/integrated-parkisons-database/

Extracted CSV files from the Integrated Parkinson’s Database should be placed in `data/raw/` to recreate the analysis

## Layout

```
data/raw/            extracted CSV files from the Integrated Parkinson’s Database
data/preprocessed/   preprocessed version of each trial, will be generated during the analysis
analysis/            numbered code files that recreate the analysis
analysis/maps/       maps that group medication, baseline characteristics and UPDRS questions
analysis/seeds/      seeds for training/validation split, will be generated during the analysis
analysis/utils/      additional functions needed for the analysis
results/             results, will be generated during the analysis
```

## Analysis

1. Create a `data/raw/` directory and place the extracted CSV files from the Integrated Parkinson’s Database in this directory
2. Run the numbered code files in `analysis/` in order of the file numbers to recreate the analysis

## Licence
Code: MIT\
Source data: governed by the data sharing agreement for the Integrated Parkinson’s Database of the Critical Path Institute, not by the MIT licence

## Funding 
This study was funded by a Strategic Basic Research Fellowship from the Research Foundation – Flanders, Belgium (FWO) for HH (FWO project number 1S30225N); by a research grant from the Research Foundation – Flanders, Belgium (FWO) for MDV and WJ (Research Project: "Artificial Intelligence (AI) for data-driven personalised medicine", FWO project number G0C9623N); by the Horizon Europe Research and Innovation Programme HORIZON-HLTH-2022-STAYHLTH: “Artificial intelligence-based Parkinson’s disease risk assessment and prognosis (AI PROGNOSIS)”, funded by the European Union, under Grant Agreement #101080581; by the Horizon Europe Research and Innovation Programme, HORIZON-HLTH-2022-IND-13: "Privacy compliant health data as a service for AI development (PHASE IV AI)”, funded by the European Union, under Grant Agreement #101095384; and received funding from the Flemish Government (AI Research Program). 