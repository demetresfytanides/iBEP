# iBEP
⸻

README — Improved BEP (Unified Multi-Layer Mixing-Length Urban Canopy Scheme)

Overview

This repository contains the implementation of the Improved BEP Urban Canopy Scheme, an enhanced version of the Building Effect Parameterization (BEP) and BEP–BEM family. The scheme introduces a unified multi-layer, vertically varying mixing-length formulation for canopy turbulence, derived from high-resolution large-eddy simulations (LES) and validated against idealized and realistic urban configurations.

The updated formulation replaces the traditional uniform mixing-length assumption within the canopy and provides a more physically realistic representation of turbulence exchange in urban street canyons. The improvements directly address long-standing limitations in predicting street-level heat exposure, humidity variability, and cooling energy demand in heterogeneous urban environments.

This package includes:
	•	WRF modifications for integrating the improved BEP scheme into the WRF Urban Physics framework.
	•	A standalone 1D column model implementation used for development, verification, and sensitivity testing.

⸻

Scientific Basis

The improved BEP scheme is grounded in:
	1.	High-resolution LES using nek5000/nekRS
LES runs covering a wide range of urban morphologies (densities, building aspect ratios, and configurations) were used to derive and verify the non-monotonic vertical structure of mixing length within the urban canopy.
	2.	Idealized and realistic urban configurations
LES was performed over canonical urban arrays.
	3.	Published data from the literature
Turbulence structure, velocity profiles, and mixing-length distributions from previous urban LES/RANS studies were integrated to ensure consistency and generality.
	4.	A unified mixing-length formulation
The resulting closure varies smoothly with height and plan area density (λₚ), ensuring:
	•	Realistic turbulence suppression in dense canopies
	•	Enhanced exchange above and near roof level
	•	Correct behavior across λₚ ∈ [0, 1]
	•	Compatibility with existing RANS and possibly RANS/LES frameworks

⸻

Contents of this Repository

iBEP/
│
├── WRF_modifications/
│   ├── <modified BEP physics routines>
│   ├── README_WRF_notes.md
│
└── 1D_implementation/
    ├── turbulence_closure.F90
    ├── example configuration files

1. WRF_modifications/

This folder contains the source-code changes needed to activate the improved BEP scheme in WRF.
All modifications are isolated to the BEP/BEP-BEM code path and retain backward compatibility.

Included:
	•	New mixing-length computation module
	•	Updated eddy diffusivity and stress calculations
	•	Optional debug flags for diagnosing canopy-layer behavior
	•	Instruction file for compilation and namelist changes

2. 1D_implementation/

This is the standalone single-column model (SCM) version used during development.
It includes:
	•	Full implementation of the unified mixing-length closure
	•	Simple force
	•	Scripts for reproducing key diagnostic figures
	•	Functions for sensitivity tests over λₚ, building height, and stability

⸻

Key Features
	•	Vertically varying mixing length derived from high-resolution LES
	•	No city-specific tuning (only morphology inputs such as λₚ and building height are required)
	•	Improved prediction of street-level temperature and humidity
	•	Better representation of vertical exchange near roof level
	•	Consistent performance across coarse (WUDAPT) and high-resolution (HiTAB) morphology data when compared against real world observations
	•	Fully compatible with the existing BEP/BEP-BEM framework
	•	Drop-in replacement for WRF urban physics modules

⸻

Usage

WRF Integration
	1.	Copy the contents of WRF_modifications/ into your WRF physics directory.
	2.	Recompile WRF.
	3.	Enable the improved BEP option via the namelist as documented in README_WRF_notes.md.

1D Model

Navigate to 1D_implementation/ compile the code using a fortran compiler.
Outputs include mixing-length profiles, turbulent fluxes, and canopy-layer diagnostics.

⸻

Citation

If you use this code, please cite:

Fytanidis, D.K. et al. (2025). How urban heterogeneity and turbulence shape street-level heat exposure. (submitted).
High-resolution LES data used for model development obtained using nek5000/nekRS.

⸻

Contact

For questions, comments, or requests for collaboration:

Dimitrios K. Fytanidis, Ph.D.
Computational Earth Scientist
Argonne National Laboratory
dfytanidis@anl.gov

⸻
