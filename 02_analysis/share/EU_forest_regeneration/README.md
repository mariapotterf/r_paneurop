# README

## Overview

This repository contains the core R scripts used for the analysis presented in the manuscript:

> **Tree regeneration after unprecedented forest disturbances in Central Europe is robust but maladapted to future climate change**  
> M. Potterf et al., 2025  
> Submitted to *Nature Ecology & Evolution*

The scripts analyze post-disturbance tree regeneration patterns across 849 field plots in 10 Central European countries and assess their climatic suitability under future climate scenarios.

All scripts are organized in a modular pipeline, from data exploration to statistical modeling and simulation outputs.


## 📁 Folder Structure

public/
├── code/ # Analysis scripts (00–05)
├── data/ # Cleaned input data
├── figs/ # Output figures
├── model/ # iLand simulation outputs
├── tables/ # Results for manuscript


---

## 📜 Script Descriptions

### `00_paths_functions.R`
Sets ups the working directory. contains the paths to teh main files. 
Utility and helper functions used across other scripts (e.g., data wrangling, plotting themes, model setup).  
➡️ Load this script at the beginning of all others.

---

### `01_analyse_structure.R`
Analyzes **field-observed regeneration structure and composition**.

**Outputs:**
- Stem density and vertical layer composition  
- Species richness and frequency across plots  
- Summary tables and plots supporting early recovery results (Figure 1, Table S1–S2)

---

### `02_model_drivers.R`
Fits **Generalized Additive Models (GAMs)** to identify climatic, soil, and disturbance predictors of stem density.

**Core analyses:**
- Tweedie-distributed GAMs (with spatial smooth)  
- Model selection via univariate AIC and drop-one tests  
- Outputs support Figures 2 and S2–S4, Table S3–S5

---

### `03_compare_adv_del_conditions_wilcox.R`
Compares **site conditions** between delayed vs. advanced regeneration using non-parametric Wilcoxon tests.

**Implements:**
- Classification into regeneration types  
- Statistical testing and visualization (Figure 3, Table S4)

---

### `04_forest_simulation.R`
Processes **iLand simulation outputs** for 30-year regeneration projections.

**Includes:**
- Infilling dynamics  
- Seed input sensitivity (Figure S5)

---

### `05_climate_suitability.R`
Evaluates **long-term climatic suitability** of regenerating species using species distribution models (SDMs).

**Calculates:**
- Species-/plot-level climatic suitability under RCP2.6, 4.5, 8.5  
- Supports Figure 4, Table 1, Table S6

---

