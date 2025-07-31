# PREFIRE_SFC

### This is a package (written primarily in MatLab, with some Python and Fortran_2003) to produce the PREFIRE 2B-SFC product. This product provides retrieved spectral surface emissivities.

This code is released under the terms of this [LICENSE](LICENSE).  The version of this package can be found in [VERSION.txt](VERSION.txt).

# Installation

## Requirements

Python version 3.6+ is required.  A recent MatLab software version (R2020b or
newer) and run-time license are also required.

The associated (Python-based) git repositories ['PREFIRE_PRD_GEN'](https://github.com/UW-PREFIRE/PREFIRE_PRD_GEN) and (Python- and Fortran-based) ['PREFIRE_pyPCRTM'](https://github.com/UW-PREFIRE) and ['PREFIRE_PCRTM_V3.4'](https://github.com/UW-PREFIRE) are also required for the proper operation of this package.

## Python Environment Setup

It is recommended to install the above Python packages in a dedicated conda environment (or something similar).  The packages used (and their versions) can be found in [conda_env.list](conda_env.list).

For example, using conda (and specifying Python 3.10.x and gfortran 8.5 from the conda-forge channel):

```
conda create --name for_PREFIRE_SFC -c conda-forge python=3.10;
conda activate for_PREFIRE_SFC;
conda install -c conda-forge gfortran=8.5;
```

The location of 'PREFIRE_pyPCRTM' depends on the value of the user's PYTHONPATH and/or sys.path -- for example, one could simply add each of those git repositories' local root Python source code directory to PYTHONPATH. Operationally, however, this package uses symbolic links to those git repositories' local root Python source code directories (or full copies of the same) in the source/ directory.

## Environment Variables

### Each job (executing this science algorithm package) is configured via information contained within environment variables.

### To specify that numpy, scipy, et cetera used by this algorithm should not use more than one thread or process, the below environment variables are expected to be set:

```
MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
```

### Some environment variables are always required to be set (also see test/run_m*.sh):

PACKAGE_TOP_DIR  :  the top-level directory (i.e., the one that contains dist/, test/, etc.) of this package

ANCILLARY_DATA_DIR  :  the package's ancillary data directory (should be an absolute path)

OUTPUT_DIR  :  the directory in which all meaningful output will be written (should be an absolute path)

ATRACK_IDX_RANGE_0BI  :  coded frame (i.e., along-track segment) subset to process and output, for example: "ATRACK_IDXRANGE_0BASED_INCLUSIVE:2001:3100" (atrack dimension indices from 2001 through 3100), "ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END" (atrack dim indices from 0 through the last frame)

PCRTM_DIR  :  the directory in which the compiled/built PCRTM (PREFIRE_PCRTM_V3.4) installation is located (contains lib/ and include/ ; should be an absolute path)

PCRTM_INPUT_DIR  :  the directory in which the ancillary PCRTM input files are located (should be an absolute path, and must end with '/')

PROC_MODE  :  If '1', create part of a 2B-SFC product granule (OE method); If '2', create part of a 2B-SFC granule (OE method + 2B-ATM prior)

L1B_RAD_FILE  :  filepath of the "source" 1B-*RAD product granule (should be an absolute path)

AUX_MET_FILE  :  filepath of the "source" AUX-MET product granule (should be an
absolute path)

L2B_MSK_FILE  :  filepath of the "source" 2B-MSK product granule (should be an
absolute path)

### Some environment variables may not need to be set for operational use (instead, some have corresponding hard-coded default values that are "baked into" each operational algorithm delivery), but exist to enable efficient development and testing (also see test/run*.sh):

PRODUCT_FULLVER  :  the full product processing/revision version string (e.g., "R01_P00").  Only increment 'Rxx' when the resulting products will be DAAC-ingested.

SRF_DISAMBIG_STR  :  identifier string that disambiguates which SRF file(s) should be used (part of the SRF filenames; e.g., "SRF_v12_2023-08-09")

# Running the test script(s)

## Obtain and unpack ancillary and test data

None (for this version).

### Prepare the test input and output directories:

`cd test;`

On Linux/UNIX systems, possibly create a useful symbolic link to the test input data (if needed):

`ln -s WHEREEVER_THE_DATA_IS/inputs inputs;`

Prepare the output directory (Linux/UNIX example):

`mkdir -p outputs;`

_OR_ perhaps something like

`ln -s /data/users/myuser/data-PREFIRE_SFC/outputs outputs;`

## Run the SFC package

### A Linux/UNIX example

`cp run.sh my-run.sh;`

Edit `my-run.sh` as needed (e.g., change input file names)

`./my-run.sh`

The output file(s) will be in `test/outputs/`

## _The creation of this code was supported by NASA, as part of the PREFIRE (Polar Radiant Energy in the Far-InfraRed Experiment) CubeSat mission._
