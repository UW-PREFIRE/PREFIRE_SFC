##!/bin/sh
#
# Intended to contain instructions to properly set up the computing environment
# (especially on large managed clusters/supercomputers).  These will be
# executed during code compilation and/or at runtime.
#
# This code has been written using "heirloom" sh-syntax, which should
#  allow the greatest (easily) possible portability among "modern" computing
#  systems.  All modifications to this code should therefore employ only
#  "heirloom" sh-syntax (e.g., *NO* syntax, assumptions, or structures
#  specific to bash, ksh, tcsh, zsh, or any nonstandard "heirloom" sh).

# Remove all currently loaded modules:
module purge;

# Load and/or unload modules:
  # If MatLab R2021b is loaded above, it requires GCC 8.x to create and use
  #  Fortran MEX executables-- so a conda environment in which gcc 8.x and
  #  gfortran 8.x are installed may be needed instead of the below module
  #  load...

