# RPDPM_MCItoDementia — Disease Progression Modelling from MCI to Dementia Using Clinical Parameters

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

This repository contains the MATLAB code used in the paper:

> **Clinical parameters predicted the progression to dementia in oldest old patients with mild cognitive impairment (MCI)**  
> Molina-Tore N., Platero C., et al.  
> *International Psychogeriatrics*, 2025.  
> DOI: [10.1016/j.inpsyc.2025.100129](https://doi.org/10.1016/j.inpsyc.2025.100129)  
> Open Access (UPM): [https://oa.upm.es/92640/](https://oa.upm.es/92640/)

---

## Overview

This toolbox applies the **Robust Parametric Disease Progression Model (RPDPM)** to characterise the long-term trajectory from **Mild Cognitive Impairment (MCI) to dementia** in oldest old patients (≥ 80 years), using exclusively **clinical parameters** collected over a **3-year prospective follow-up** from a proprietary patient cohort.

The RPDPM framework linearly maps individual ages to a **Disease Progression Score (DPS)** and jointly fits constrained generalised logistic functions to the longitudinal dynamics of biomarkers using M-estimation, providing robustness against outliers and missing data. The estimated inflection points are used to temporally order clinical markers along the disease course, and the resulting DPS is used for clinical status classification via a Bayesian classifier.

Key features of this implementation:
- Modelling of MCI-to-dementia progression using **clinical parameters only** (no neuroimaging required)
- Application to **oldest old patients** (≥ 80 years), an underrepresented population in DPM literature
- Longitudinal cohort with **3-year follow-up**
- Multi-start strategy (`multi_RPDPM_MCItoDementia.m`) for robust parameter estimation
- Outputs include fitted biomarker trajectories, DPS estimates, and classification results

---

## Repository Structure

```
RPDPM_MCItoDementia/
├── multi_RPDPM_MCItoDementia.m   # Main script: multi-start RPDPM fitting and evaluation
├── rpdpm/                        # Core RPDPM algorithm functions
├── aux_dpm/                      # Auxiliary functions (preprocessing, plotting, evaluation)
├── data/                         # Placeholder for input data (not distributed — see Data section)
├── output/                       # Results: figures, tables, fitted parameters
└── README.md
```

---

## Requirements

- **MATLAB** ≥ R2020a
- **Statistics and Machine Learning Toolbox**
- **Optimization Toolbox**

No additional external toolboxes are required.

---

## Data

The experiments in this paper use data from a **proprietary longitudinal cohort** of oldest old patients (≥ 80 years) with MCI, followed over 3 years at a clinical centre. Due to patient privacy and ethical constraints, the data cannot be redistributed as part of this repository.

The `data/` folder contains a **synthetic example dataset** with the same format as the original data, allowing users to run the pipeline end-to-end and verify the implementation.

Clinical parameters used as biomarkers include cognitive assessments (MMSE, MoCA, CDR-SB), functional scales, and other routinely collected clinical measures. See the paper for full details.

---

## Usage

Open MATLAB, navigate to the root of the repository, and run the main script:

```matlab
% From the MATLAB command window or editor:
cd /path/to/RPDPM_MCItoDementia
multi_RPDPM_MCItoDementia
```

The script will:
1. Load input data from `data/`
2. Run the multi-start RPDPM fitting procedure
3. Estimate individual Disease Progression Scores (DPS)
4. Fit generalised logistic trajectories to each clinical biomarker
5. Classify clinical status using a Bayesian classifier
6. Save figures and result tables to `output/`

Configurable parameters (number of random initialisations, biomarkers to model, train/test split) are set at the top of `multi_RPDPM_MCItoDementia.m`.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{molinatore2025clinical,
  title   = {Clinical parameters predicted the progression to dementia in oldest old patients
             with mild cognitive impairment ({MCI})},
  author  = {Molina-Tore, N. and Platero, Carlos and others},
  journal = {International Psychogeriatrics},
  year    = {2025},
  doi     = {10.1016/j.inpsyc.2025.100129},
  url     = {https://doi.org/10.1016/j.inpsyc.2025.100129}
}
```

The RPDPM algorithm implemented here was originally proposed in:

```bibtex
@article{ghazi2021robust,
  title   = {Robust parametric modeling of {Alzheimer's} disease progression},
  author  = {Ghazi, Mostafa Mehdipour and Nielsen, Mads and Pai, Akshay and
             Cardoso, M. Jorge and Leung, Kelvin K. and Ourselin, S{\'e}bastien
             and Sommer, Stefan},
  journal = {NeuroImage},
  volume  = {225},
  pages   = {117460},
  year    = {2021},
  doi     = {10.1016/j.neuroimage.2020.117460}
}
```

---

## Related Repositories

Other DPM-related tools from the same research group:

- [preAD_DPM](https://github.com/cplatero/preAD_DPM) — Benchmarking parametric DPMs (GRACE, Leaspy, RPDPM) for early AD detection
- [twogrsurvana](https://www.nitrc.org/projects/twogrsurvana/) — Longitudinal survival analysis and two-group comparison for MCI-to-AD progression
- [predict_mci2ad](https://www.nitrc.org/projects/predict_mci2ad/) — Predicting Alzheimer's conversion in MCI patients using neuroimaging and clinical markers

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

Carlos Platero — Universidad Politécnica de Madrid  
[https://github.com/cplatero](https://github.com/cplatero)
