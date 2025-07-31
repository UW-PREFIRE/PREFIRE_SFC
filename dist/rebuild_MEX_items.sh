#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./rebuild_MEX_items.sh   OR   bash rebuild_MEX_items.sh

# Driver for the recompile/rebuild of all MEX (Fortran-Matlab interface with
#  PCRTM) items for this science algorithm package.

absfpath() {
  # Generate absolute filepath from a relative (or even an absolute) filepath.
  #
  # Based on (circa Oct 2023) https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
  # 
  # $1     : a relative (or even an absolute) filepath
  # Returns the corresponding absolute filepath.
  if [ -d "$1" ]; then
    # dir
    (cd "$1"; pwd)
  elif [ -f "$1" ]; then
    # file
    if [[ $1 = /* ]]; then
      echo "$1"
    elif [[ $1 == */* ]]; then
      echo "$(cd "${1%/*}"; pwd)/${1##*/}"
    else
      echo "$(pwd)/$1"
    fi
  fi
}

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package dist/ directory)
readonly base_dir="$(absfpath ".")";

this_top_dir="$(absfpath "${base_dir}/..")";

# Execute any necessary machine setup instructions ((un)loading modules, etc.):
. "${this_top_dir}/dist/perform_machine_setup.sh";

# Change dir to the 'dist' directory:
begdir=`pwd`;  # Save initial dir
cd "${this_top_dir}/dist";

# Set input paths within environment variables:
PCRTM_DIR="${this_top_dir}/dist/subpackages/PREFIRE_PCRTM_V3.4";
#PCRTM_DIR=/data/rttools/PCRTM/gfortran_8.5_build/PREFIRE_PCRTM_V3.4;

THIS_F_SRCDIR="${this_top_dir}/source/fortran";
BUILD_DIR="${this_top_dir}/source/matlab";

GFORTRAN_PARADIGM='8.x'
#GFORTRAN_PARADIGM='11.x'

export PCRTM_DIR THIS_F_SRCDIR GFORTRAN_PARADIGM BUILD_DIR;

# Remove old MEX-related detritus:
rm -f *.mod *.mexa64;

# Build MEX items:
matlab -singleCompThread -nodisplay -batch app_MEX_build;

# Return to initial dir:
cd "${begdir}";

# Check for existence of PCRTM_MEX_interface.mexa64
tmpf="${BUILD_DIR}/PCRTM_MEX_interface.mexa64";
test -f "${tmpf}" || { echo "ERROR: ${tmpf} was not created."; exit 1; }

echo 'Finished (re)building MEX items.';
