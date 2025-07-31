#include "fintrf.h"

! This is a MEX file for MATLAB.

!------------------------------------------------------------------------------
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!+ Gateway routine.

use SFC_driver_for_PCRTM_mod, only: SFC_driver_for_PCRTM

implicit none

! mexFunction arguments:
mwPointer, dimension(*) :: plhs, prhs
integer :: nlhs, nrhs

! Function declarations:
mwPointer :: mxCreateDoubleMatrix, mxGetString, mxIsInt8

#if MX_HAS_INTERLEAVED_COMPLEX
mwPointer :: mxGetDoubles, mxGetInt8s
#else
mwPointer :: mxGetPr
#endif

integer :: mxIsNumeric, mxIsChar
mwPointer :: mxGetM, mxGetN

! Pointers to input/output mxArrays:
mwPointer :: x_pr, ans1_pr, ans2_pr

! Array information:
mwPointer :: m, n
mwSize :: size_of_this

! Declare arguments for computational routine:
integer, parameter :: numel=20000  ! Max array size
real(8), dimension(numel) :: f8_flat_input, atm_1Dfields, cld_1Dfields, emis
integer(1), dimension(numel) :: i1_flat_input
character(len=120) :: input_cbuffer
real(8), dimension(numel) :: ans1, ans2

integer :: i, n_atm_1dfields, n_cld, n_cld_1Dfields, n_e1, n_e2, n_f8fi,   &
     n_i1fi, n_lev, status, ic, n_a1

! Check for proper number of arguments. 
if (nrhs /= 6) then
   call mexErrMsgIdAndTxt('MATLAB:SFC_driver:nInput',  &
        'Six inputs required (i1_flat_input, f8_flat_input, atm_1Dfields, cld_1Dfields, emis, PCRTM_static_input_dir).')
else if (nlhs /= 2) then
   call mexErrMsgIdAndTxt('MATLAB:SFC_driver:nOutput', 'Two outputs required.')
end if

! Get the size(s) of the input array(s):
do i=1,nrhs
   m = mxGetM(prhs(i))
   n = mxGetN(prhs(i))
   size_of_this = m*n

   if (size_of_this > numel) then
      call mexErrMsgIdAndTxt ('MATLAB:SFC_driver:mSize',  &
           'row * column must be <= 1000')
   end if

   if (i <= 5) then
      ! Check that the array is numeric (not strings):
      if (mxIsNumeric(prhs(i)) == 0) then
         print *,'NonNumeric ERROR: ', i
         call mexErrMsgIdAndTxt('MATLAB:SFC_driver:NonNumeric',  &
              'Input must be a numeric array.')
      end if
   else
      ! Check that this is string data:
      if (mxIsChar(prhs(i)) /= 1) then
         print *,'NonString ERROR: ', i
         call mexErrMsgIdAndTxt('MATLAB:SFC_driver:NonString',  &
              'Input must be a string.')
      end if
   end if

   if (i >= 2 .and. i <=5) then  ! 8-byte floating-point data
#if MX_HAS_INTERLEAVED_COMPLEX
      x_pr = mxGetDoubles(prhs(i))
#else
      x_pr = mxGetPr(prhs(i))
#endif
      ! Create a Fortran array from the input, which is the mxArray 'x_pr':
      select case(i)
      case(2)
         n_f8fi = size_of_this
         call mxCopyPtrToReal8(x_pr, f8_flat_input, size_of_this)
      case(3)
         n_lev = m; n_atm_1Dfields = n
         call mxCopyPtrToReal8(x_pr, atm_1Dfields, size_of_this)
      case(4)
         n_cld = m; n_cld_1Dfields = n
         call mxCopyPtrToReal8(x_pr, cld_1Dfields, size_of_this)
      case(5)
         n_e1 = m; n_e2 = n
         call mxCopyPtrToReal8(x_pr, emis, size_of_this)
      end select

   else if (i == 1) then  ! 1-byte integer data
      n_i1fi = size_of_this
      ! Check for proper input type:
      if (mxIsInt8(prhs(i)) == 0) then
         call mexErrMsgIdAndTxt('MATLAB:SFC_driver:WrongType',  &
              'Input must be int8 array.')
      end if
#if MX_HAS_INTERLEAVED_COMPLEX
      x_pr = mxGetInt8s(prhs(i))
#else
      x_pr = mxGetPr(prhs(i))
#endif
      ! Create a Fortran array from the input, which is the mxArray 'x_pr':
      call mxCopyPtrToInteger1(x_pr, i1_flat_input, size_of_this)

   else  ! Character data (PCRTM_static_input_dir)
      ! The input must be a row vector:
      if (m /= 1) then
         call mexErrMsgIdAndTxt('MATLAB:SFC_driver:NonRowVector',  &
              'Input must be a row vector.')
      end if
      status = mxGetString(prhs(i), input_cbuffer, 120)
        ! Add required separator, if not present:
      ic = len_trim(input_cbuffer)
      if (input_cbuffer(ic:ic) /= '/') input_cbuffer(ic+1:ic+1) = '/'
   end if
end do

!call mexMakeMemoryPersistent(ptr)

! Call the computational subroutine:
call SFC_driver_for_PCRTM(n_i1fi, i1_flat_input, n_f8fi, f8_flat_input,  &
                          n_lev, n_atm_1Dfields, atm_1Dfields,  &
                          n_cld, n_cld_1Dfields, cld_1Dfields,  &
                          n_e1, n_e2, emis, input_cbuffer,  &
                          numel, n_a1, ans1, ans2)

! Create matrix for the return arguments:
plhs(1) = mxCreateDoubleMatrix(n_a1, 1, 0)
plhs(2) = mxCreateDoubleMatrix(n_a1, 1, 0)

#if MX_HAS_INTERLEAVED_COMPLEX
   ans1_pr = mxGetDoubles(plhs(1))
   ans2_pr = mxGetDoubles(plhs(2))
#else
   ans1_pr = mxGetPr(plhs(1))
   ans2_pr = mxGetPr(plhs(2))
#endif

! Load the data into ans*_pr, which is the output to MATLAB.
call mxCopyReal8ToPtr(ans1, ans1_pr, n_a1)
call mxCopyReal8ToPtr(ans2, ans2_pr, n_a1)

end subroutine mexFunction
