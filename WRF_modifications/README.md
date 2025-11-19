# Improved BEP (iBEP) — Unified Multi-Layer Mixing-Length Urban Canopy Scheme
### WRF drop-in replacement module

This repository contains a single modified WRF physics file:

- **module_sf_bep_bem.f90**  
  An improved version of the BEP/BEP–BEM urban canopy scheme implementing a unified, vertically varying mixing-length formulation for urban canopy turbulence.  
  The closure is derived from high-resolution LES (nek5000/nekRS) and validated against both idealized and realistic urban configurations.

The file is a direct drop-in replacement for the standard BEP–BEM module shipped with WRF.

## Usage

1. Copy the provided file into your WRF physics directory:
   ```bash
   cp module_sf_bep_bem.f90 /path/to/WRF/phys/module_sf_bep_bem.f90
   ```
   Replace the existing version.

2. Compile WRF as you normally do.  
   No special steps are required.

## Namelist Configuration

To activate the scheme in WRF, enable BEP–BEM in namelist.input:
```fortran
sf_urban_physics = 2
```

No additional flags or custom settings are required beyond your normal urban configuration.

## Model Notes

- iBEP uses the standard morphology inputs used by BEP/BEM (plan-area fraction λₚ, building height, roof/wall fractions).
- No city-specific parameters or tuning are included.
- The scheme improves turbulence representation inside the canopy, vertical exchange near roof level, and resulting temperature/RH patterns.

## Citation

If you use this module, please cite:

Fytanidis, D.K. et al. (2025). *How urban heterogeneity and turbulence shape street-level heat exposure.* (Submitted).  
High-resolution LES development supported by nek5000/nekRS.

## Contact

For questions or collaboration:

Dimitrios K. Fytanidis, Ph.D.  
Computational Earth Scientist  
Argonne National Laboratory  
dfytanidis@anl.gov
