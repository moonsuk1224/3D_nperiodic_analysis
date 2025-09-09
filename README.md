# 3D N-Periodic SSFP Analysis

This repository contains code for **3D reconstruction and analysis of N-periodic Steady-State Free Precession (SSFP) MRI data**, developed as part of my MSc project at Imperial College London.  
It is adapted and extended from my academic supervisor [Pete Lally’s `nperiodic_analysis`](https://github.com/petelally/nperiodic_analysis.git).

---

## 🔹 Overview
The project aims to:
- Implement **3D N-periodic SSFP reconstruction**, extending the original 2D framework.
- Enable comparison with **multi-echo FLASH (meFLASH)** acquisitions for **Quantitative Susceptibility Mapping (QSM)**.
- Investigate the feasibility of SSFP as an alternative to GRE for brain QSM by leveraging its high SNR efficiency and multi-echo properties.

---

## 🔹 Features
- **Raw Siemens `.dat` file support** (via [`mapVBVD`](https://github.com/MarkusDietrich/mapVBVD)).
- **3D POCS partial Fourier reconstruction** for N-periodic SSFP.
- **F-state demodulation** and separation from phase-cycled data.
- **Coil combination** via SOS, ESPIRiT, and SENSE.
- Export to **DICOM/NIfTI** for downstream QSM analysis (e.g., [SEPIA toolbox](https://github.com/kschan0214/sepia)).

---

## 🔹 Repository Structure
/scripts → main reconstruction scripts
/external → external functions (e.g. mapVBVD, ESPIRiT)
/figures → example figures and schematic diagrams
/data_dummy → placeholder folder (real MRI data not included)
/utils → helper functions (DICOM export, k-space cropping, etc.)

