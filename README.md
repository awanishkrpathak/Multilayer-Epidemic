# Computational Codes

This repository contains the Fortran programs and Jupyter notebooks used for the
computational analyses.

## Repository structure

### `Fortran_Codes/`
- `influenza_lyon_final.f90`
- `RSV_Lyon_final.f90`
- `SarCov_lyon_final.f90`
- `Mozambique_final.f90`

### `Notebooks/`
- `Lyon_Hospital/Lyon_Hospital_NatCom.ipynb`
- `Mozambique/Mozambique_NatCom.ipynb`
- `Additional_Analysis/Natcom_Final.ipynb`

The notebooks include the Python analysis and plotting code associated with the
corresponding computational workflows.

## Running the notebooks

Install Python and JupyterLab, then install the scientific Python packages used
by the notebooks. Open the required notebook in JupyterLab and run its cells in
order. Any input files referenced by a notebook must be available at the paths
expected by that notebook.

## Running the Fortran programs

A Fortran compiler such as GNU Fortran (`gfortran`) is required. Compile and run
the relevant source file using the compiler options and input files required by
that program. The precise execution commands and input requirements should be
checked against each source file before running.

## Reproducibility

The files in this repository are provided as the computational source materials
for the study. The notebooks have been checked for valid notebook structure and
Python syntax, but the full numerical workflows have not been executed as part
of packaging. Please consult the associated manuscript for the model description
and parameter definitions.
