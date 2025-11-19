# 1D iBEP Column Model

This directory contains a 1D single-column implementation of the improved BEP (iBEP) mixing-length urban canopy scheme. The code solves for vertical profiles of velocity, turbulence, and scalars in an idealized urban canopy column using the unified multi-layer mixing-length formulation.

Implementation based on code by Dr A. Martilli see for example:
Santiago, J. L., & Martilli, A. (2010). A dynamic urban canopy parameterization for mesoscale models based on computational fluid dynamics Reynolds-averaged Navier–Stokes microscale simulations. Boundary-layer meteorology, 137(3), 417-439. 

## Files
- `iBEP.f90` – Fortran source code for the 1D column model.
- `input_column` – Input file defining the vertical grid, forcing, and morphology parameters.

## Requirements
- A Fortran compiler (e.g., `gfortran`).
- On macOS ARM (M1/M2/M3), use an ARM-native compiler, for example:
  ```bash
  conda install -c conda-forge gfortran_osx-arm64
  ```

## Compilation
From inside this directory:
```bash
gfortran -O3 -o iBEP iBEP.f90
```
```

## Running the Model
Ensure `input_column` is in the same directory, then run:
```bash
./iBEP
```

The model automatically reads `input_column`.

## Output
The model writes:
- `output` – main diagnostic file
- Additional files such as:
  - `iBEPoutput_staggered_r`
  - `iBEPoutput_aligned_r`

depending on the configuration.

## Contact
Dimitrios K. Fytanidis  
dfytanidis@anl.gov
