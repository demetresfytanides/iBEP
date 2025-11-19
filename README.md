# iBEP — a Unified Multi-Layer Mixing-Length Urban Canopy Scheme

## Overview

This repository contains the implementation of the Unified Multi-Layer Mixing-Length Urban Canopy Scheme (iBEP), an enhanced version of the standard Building Effect Parameterization (BEP) and BEP–BEM family.

iBEP introduces a unified multi-layer, vertically varying mixing-length formulation for canopy turbulence, derived from high-resolution large-eddy simulations (LES) using nek5000/nekRS and validated against both idealized and realistic urban configurations.

The updated formulation replaces the traditional uniform mixing-length assumption within the canopy, providing a more physically realistic representation of turbulent exchange in urban street canyons. This improvement directly addresses long-standing limitations in predicting street-level heat exposure, humidity variability, and cooling energy demand in heterogeneous urban environments.

This repository includes:
- WRF modifications implementing the new BEP module.
- A standalone 1D column model used for development, verification, and sensitivity testing.

---

## Scientific Basis

The new BEP scheme is grounded in:

1. **High-resolution LES using nek5000/nekRS**  
   LES simulations spanning a wide range of urban morphologies—plan-area density, aspect ratio, and building configuration—were used to derive and verify the non-monotonic vertical structure of mixing length within the canopy.

2. **Idealized and realistic urban configurations**  
   The formulation was evaluated against canonical cube arrays (aligned and staggered) and realistic neighborhood configurations.

3. **Published literature and datasets**  
   Turbulence structure, velocity profiles, and mixing-length distributions from earlier LES/RANS studies were incorporated to ensure consistency and generality.

4. **Unified mixing-length formulation**  
   The closure:
   - varies smoothly with height,
   - depends on plan-area density (λₚ),
   - reproduces turbulence suppression in dense canopies,
   - enhances vertical exchange near and above roof level,
   - behaves correctly across the full λₚ range tested,
   - and remains fully compatible with RANS and hybrid RANS/LES canopy representations.

---

## Repository Structure

	•	WRF_modifications/
	   •	Modified BEP urban physics routines
	   •	README.md (WRF-specific usage and namelist notes)
	•	1D_implementation/
	   •	Turbulence closure routines
	   •	Example configuration files
	
---

## WRF_modifications/

Contains all modifications necessary to activate the new BEP scheme in WRF.  
Changes are constrained to the BEP/BEP-BEM code path and maintain backward compatibility.

Included:
- New mixing-length computation module.
- Updated eddy diffusivity and stress calculations.
- Debug flags for diagnosing canopy-layer behavior.
- Instructions for compilation and namelist configuration.

---

## 1D_implementation/

Standalone single-column model (SCM) version of iBEP used for testing and verification.

Features:
- Complete implementation of the unified mixing-length closure.
- Forcing options for sensitivity exploration.
- Diagnostic tools for mixing-length profiles, turbulent fluxes, and canopy-layer structure.

---

## Key Features

- Vertically varying mixing length derived from high-resolution LES.
- No city-specific tuning (requires only local morphology inputs such as λₚ and building height).
- Improved prediction of street-level temperature and humidity.
- Better representation of vertical exchange near roof level.
- Strong performance for both WUDAPT and HiTAB morphology datasets.
- Fully compatible with the BEP/BEP-BEM framework.
- Drop-in replacement for existing WRF urban physics modules.

---

## Usage

### WRF Integration

1. Copy the contents of `WRF_modifications/` into your WRF physics directory.
2. Recompile WRF.
3. Enable iBEP in the WRF namelist (see `README_WRF_notes.md`).

### 1D Model

1. Navigate to `1D_implementation/`.
2. Compile using a Fortran compiler.
3. Run example scripts to generate diagnostic profiles and plots.

---

## Citation

If you use this code, please cite:

Fytanidis, D. K., Tan, H., Kotamarthi, R., Wang, J., Martilli, A., O’Brien, J., Muradyan, P., Collis, S., & Negri, C. (2025).  *How urban heterogeneity and turbulence shape street-level heat exposure.* (Manuscript submitted for publication).

---

## Contact

For questions, comments, or collaboration requests:

**Dimitrios K. Fytanidis, Ph.D.**  
Computational Earth Scientist  
Argonne National Laboratory  
dfytanidis@anl.gov
