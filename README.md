# EffiMultiPSP

This repository contains the ** R code implementations** accompanying the paper:

> **"Efficiency of multivariate tests in trials in progressive supranuclear palsy"**  
> Published in *Scientific Reports*, Nature Portfolio, 26 October 2024.  
> [Link to article](https://www.nature.com/articles/s41598-024-76668-4)

---

## Overview

The work in this repository is based on the published study, which investigates statistical approaches to analyze treatment effects in progressive supranuclear palsy (PSP) clinical trials based on a modification of the Progressive Supranuclear Palsy Rating Scale containing 10 items (PSPRS-10).  

The study compared:
- **Classical approaches**, such as PSPRS-10 sum scores.  
- **Multivariate tests** and **multiple comparisons approaches** for item-level analyses.  
- **Novel Item Response Theory (IRT)–based models** to measure disease status.  

Extensive simulation studies were performed to evaluate these approaches, and a re-analysis of the ABBV-8E12 clinical trial was conducted to illustrate their use.  

**Key findings:**  
- PSPRS sum scores are robust but their efficiency depends on the treatment effect pattern.  
- IRT-based methods provide the highest power when data are generated from IRT models.  
- Multiple testing–based approaches are advantageous when treatment effects are localized to specific items/domains.  
- No single test is universally optimal — the choice depends on the effect size patterns.  

---

## Repository Structure

- **`R_MainManuscript/`**  
  Contains the main R implementations used for the manuscript simulations and analyses.  
  ⚠️ **Note:** These scripts were implemented using **real clinical trial data**, which **cannot be made publicly available**.  

- **`effMultiend/`**  
  An R package under development.  
  - Currently focuses on **simulation studies using an Item Response Theory (IRT) model** (the third data-generating mechanism described in the paper).  
  - Future updates will extend the package to include additional simulation methods and analysis approaches.  

- **Other supporting files**  
  `.gitignore`, `.Rbuildignore`, `DESCRIPTION`, `NAMESPACE`, and project metadata for package development.

- **Note:** This repository does not yet include the additional simulation codes described in *Supplementary Material A* of the paper. These materials are currently private but will be made publicly available and added here in a future update.

---

## Status

- The `R_MainManuscript` folder provides the raw code basis for the published paper.  
- The `effMultiend` package is **under development** and not yet intended for release.

---

## Citation

If you use this code or package in your research, please cite:

Yousefi, E., et al. (2024). *Efficiency of multivariate tests in trials in progressive supranuclear palsy*. Scientific Reports.  
[https://doi.org/10.1038/s41598-024-76668-4](https://doi.org/10.1038/s41598-024-76668-4)

---

## License

This repository is released under the GNU License. See the [LICENSE](LICENSE) file for details.
