---

## WRF_modifications/

Contains all modifications necessary to activate the improved BEP scheme in WRF.  
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

Fytanidis, D.K. et al. (2025). *How urban heterogeneity and turbulence shape street-level heat exposure.* (Submitted).  
High-resolution LES data were obtained using nek5000/nekRS.

---

## Contact

For questions, comments, or collaboration requests:

**Dimitrios K. Fytanidis, Ph.D.**  
Computational Earth Scientist  
Argonne National Laboratory  
dfytanidis@anl.gov
