# UPF and Sleep Quality and Multi-Omics Analysis Project

This project contains analysis code for examining the association between ultra-processed food (UPF) consumption and sleep quality, as well as multi-omics analysis workflows.

## Study Design Overview

This study is based on two large prospective cohorts (UK Biobank and HRS) to systematically evaluate the long-term association between UPF consumption and sleep quality. By integrating proteomic and metabolomic data, we identified molecular signatures associated with UPF intake. We further analyzed the mediating roles of these molecular signatures between UPF and sleep quality, as well as their associations with multiple chronic diseases and mortality.

## File Descriptions

1. `upf_sleep_association.R` - UPF and sleep quality association analysis
   - Association analysis between UPF intake and sleep quality in the UK Biobank cohort
   - Analysis of PSQI components
   - Validation analysis in the HRS cohort

2. `omics_signature_analysis.R` - Multi-omics signature construction and analysis
   - Univariate analysis to screen molecules significantly associated with UPF
   - Elastic Net analysis to construct proteomic and metabolomic signatures
   - Functional enrichment analysis
   - Association analysis between multi-omics signatures and sleep quality

3. `outcome_wide_analysis.R` - Outcome-wide analysis
   - Association analysis between UPF intake and related multi-omics signatures with multiple disease outcomes