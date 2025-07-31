% (Re-)compile/build all MEX (Fortran-Matlab interface with PCRTM) items for
%  this science algorithm package.

% Get needed information from environment variables:
PCRTM_dir = getenv('PCRTM_DIR');
this_F_srcdir = getenv('THIS_F_SRCDIR');
gfortran_paradigm = getenv('GFORTRAN_PARADIGM');
build_dir = getenv('BUILD_DIR');

cd(build_dir);

% Uses a set of special CPU-specific options that is common between system the
%  algorithm delivery package is created on and the SDPS:
if strcmp(gfortran_paradigm, '11.x')
   fflags_str = 'FFLAGS=''$FFLAGS -fallow-invalid-boz -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -mavx -mavx2 -mf16c -mfma -mbmi -mbmi2 -march=x86-64''';
else  % '8.x'
   fflags_str = 'FFLAGS=''$FFLAGS -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -mavx -mavx2 -mf16c -mfma -mbmi -mbmi2 -march=x86-64''';
end

PCRTM_incdir = [PCRTM_dir, '/include']
PCRTM_libdir = [PCRTM_dir, '/lib']

% To "bootstrap" the process (creates a needed *.mod file, but produces an
%  incorrectly-named MEX executable):
eval(['mex -v ', fflags_str, ' -I',PCRTM_incdir, ' ', ...
      this_F_srcdir,'/SFC_driver_for_PCRTM.f90 ', ...
      this_F_srcdir,'/PCRTM_MEX_interface.F90 ', PCRTM_libdir,'/libPCRTM.a'])

% To produce the correctly-named MEX executable:
eval(['mex -v ', fflags_str, ' -I',PCRTM_incdir, ' ', ...
      this_F_srcdir,'/PCRTM_MEX_interface.F90 ',  ...
      this_F_srcdir,'/SFC_driver_for_PCRTM.f90 ', PCRTM_libdir,'/libPCRTM.a'])
