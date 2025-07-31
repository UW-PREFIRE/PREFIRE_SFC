module SFC_driver_for_PCRTM_mod

use PCRTM_CONSTANTS, only: SINGLE, DOUBLE, pbnd_101
use PCRTM_ATM_ABSORPTION_DEFINE, only: PCRTM_ATM_ABSORPTION_TYPE,  &
     PCRTM_ATM_ABS_STRUCT_TYPE
use PCRTM_TR_SOLUTION, only: PCRTM_TR_SOLUTION_TYPE
use PCRTM_CLOUD_LUT_IO, only: PCRTM_CLD_TABLE_DEF
use PCRTM_PC_SOLUTION, only: PCRTM_EOF_SOLUTION_TYPE, TRUNCATE_PCRTM_EOF_JACOB
use PCRTM_RT_SOLUTION_DEFINE, only: PCRTM_RT_SOLUTION_TYPE
use PCRTM_JACOBIAN, only: PCRTM_NM_JACOBIAN_TYPE, PCRTM_CH_JACOBIAN_TYPE,  &
     PCRTM_PC_JACOBIAN_TYPE
use PCRTM_ATMOSPHERE_DEFINE, only: PCRTM_ATMOSPHERE_TYPE
use PCRTM_ATMOSPHERE_LAYER, only: PCRTM_GEOMETRY_TYPE
use PCRTM_FORWARD_MODEL, only: PCRTM_FORWARD_RT

implicit none

type(PCRTM_ATM_ABSORPTION_TYPE) :: atm_abs_coef, atm_abs_coef_Ja
type(PCRTM_TR_SOLUTION_TYPE), allocatable, dimension(:) :: tr_solution,  &
     tr_solution_Ja
type(PCRTM_CLD_TABLE_DEF) :: ice_grid, wat_grid, ice_grid_Ja, wat_grid_Ja
type(PCRTM_EOF_SOLUTION_TYPE), allocatable, dimension(:) :: EOF_solution,  &
     EOF_solution_Ja
type(PCRTM_RT_SOLUTION_TYPE) :: RT_solution, RT_solution_Ja
type(PCRTM_NM_JACOBIAN_TYPE) :: K_M, K_M_Ja
type(PCRTM_CH_JACOBIAN_TYPE), allocatable, dimension(:) :: K_CH, K_CH_Ja
type(PCRTM_PC_JACOBIAN_TYPE), allocatable, dimension(:) :: K_PC, K_PC_Ja

type(PCRTM_ATMOSPHERE_TYPE) :: atm, atm_Ja
type(PCRTM_GEOMETRY_TYPE) :: geometry

type(PCRTM_ATM_ABS_STRUCT_TYPE) :: PCRTM_STND, PCRTM_STND_Ja

integer, parameter :: i4 = selected_int_kind(9),  &    ! can hold +/- 10^9
                      i8 = selected_int_kind(18),  &   ! can hold +/- 10^18
                      b4 = 4
integer(i4), parameter :: f4 = SINGLE, f8 = DOUBLE  ! More-convenient names

integer(i4), allocatable, dimension(:) :: nPCbnd, nPCbnd_Ja

integer(i4), parameter :: max_n_abs_molec = 6

integer(i4), parameter ::  &
          sensor_ID = 2  ! sensor identifier (PCRTM ID for CLARREO 0.5)

real(f4), dimension(size(pbnd_101)) :: CO2_0 =  &  ! Default CO2 profile [ppm]
     (/ 363.5187, 363.5187, 363.5187, 363.5187, 363.5187,  &
        363.5312, 363.5702, 363.6875, 363.8089, 363.9016,  &
        363.9719, 364.0410, 364.0981, 364.1384, 364.1780,  &
        364.2169, 364.3329, 364.5038, 364.6715, 364.8359,  &
        364.9962, 365.1496, 365.2998, 365.4469, 365.5909,  &
        365.6422, 365.6590, 365.6754, 365.6914, 365.7214,  &
        365.7626, 365.8028, 365.8419, 365.9814, 366.1763,  &
        366.3661, 366.5508, 366.7858, 367.0166, 367.2409,  &
        367.4904, 367.7629, 368.0271, 368.3109, 368.6626,  &
        369.0027, 369.3418, 369.7074, 370.0602, 370.4013,  &
        370.7320, 371.0502, 371.3610, 371.6623, 371.9639,  &
        372.4423, 372.9003, 373.5883, 374.3813, 375.0573,  &
        375.5598, 376.0378, 376.4460, 376.6990, 376.9385,  &
        377.1650, 377.3794, 377.5828, 377.7740, 377.9535,  &
        378.1273, 378.2934, 378.4482, 378.5996, 378.7636,  &
        378.9151, 379.0546, 379.1853, 379.3575, 379.5143,  &
        379.6564, 379.7847, 379.8999, 380.1389, 380.3916,  &
        380.6143, 380.8089, 381.0880, 381.6050, 382.0438,  &
        382.4115, 383.0423, 383.8282, 384.4527, 384.9500,  &
        385.3369, 385.6151, 385.7570, 385.8480, 385.8480,  &
        385.8480 /)

real(f4), dimension(size(pbnd_101)) :: N2O_0 =  &  ! Default N2O profile [ppm]
     (/ 1.0940E-04, 1.4211E-04, 2.0891E-04, 3.0890E-04, 4.0362E-04,  &
        4.9747E-04, 5.4515E-04, 5.7810E-04, 5.9339E-04, 6.3985E-04,  &
        8.3545E-04, 1.1641E-03, 1.5052E-03, 1.7441E-03, 1.8774E-03,  &
        1.9921E-03, 2.0930E-03, 2.2095E-03, 2.3661E-03, 2.5988E-03,  &
        2.8556E-03, 3.1843E-03, 3.7487E-03, 5.0795E-03, 7.0230E-03,  &
        1.1365E-02, 1.7546E-02, 2.6933E-02, 3.8778E-02, 5.2639E-02,  &
        6.8569E-02, 8.5894E-02, 1.0431E-01, 1.2393E-01, 1.4152E-01,  &
        1.6010E-01, 1.7497E-01, 1.9052E-01, 2.0179E-01, 2.1305E-01,  &
        2.2369E-01, 2.3448E-01, 2.4457E-01, 2.5434E-01, 2.6319E-01,  &
        2.6992E-01, 2.7695E-01, 2.8183E-01, 2.8675E-01, 2.9113E-01,  &
        2.9444E-01, 2.9787E-01, 3.0032E-01, 3.0230E-01, 3.0434E-01,  &
        3.0573E-01, 3.0697E-01, 3.0825E-01, 3.0921E-01, 2.9705E-01,  &
        3.0864E-01, 3.1982E-01, 3.1482E-01, 3.0030E-01, 3.1977E-01,  &
        3.2718E-01, 3.2113E-01, 3.2503E-01, 3.2076E-01, 3.1124E-01,  &
        3.1376E-01, 3.1814E-01, 3.2711E-01, 3.0301E-01, 3.0637E-01,  &
        3.0329E-01, 3.1982E-01, 3.2701E-01, 3.1335E-01, 3.2132E-01,  &
        3.1009E-01, 3.1746E-01, 3.1435E-01, 3.1964E-01, 3.0253E-01,  &
        3.0242E-01, 3.0021E-01, 3.1795E-01, 3.2329E-01, 3.0075E-01,  &
        3.0867E-01, 3.0852E-01, 3.0830E-01, 3.2395E-01, 3.1714E-01,  &
        3.2099E-01, 3.2697E-01, 3.0389E-01, 3.1763E-01, 3.1489E-01,  &
        3.0853E-01 /)

real(f4), dimension(size(pbnd_101)) :: CO_0 =  &  ! Default CO profile [ppm]
     (/ 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02, 4.2239E-02,  &
        4.2239E-02, 4.2239E-02, 4.2253E-02, 4.2291E-02, 4.2366E-02,  &
        4.2492E-02, 4.2691E-02, 4.2987E-02, 4.3407E-02, 4.3978E-02,  &
        4.4738E-02, 4.5713E-02, 4.6933E-02, 4.8405E-02, 5.0122E-02,  &
        5.2066E-02, 5.4222E-02, 5.6581E-02, 5.9133E-02, 6.1829E-02,  &
        6.4596E-02, 6.7360E-02, 7.0086E-02, 7.2749E-02, 7.5311E-02,  &
        7.7702E-02, 7.9842E-02, 8.1682E-02, 8.3178E-02, 8.4280E-02,  &
        8.4970E-02, 8.5232E-02, 8.5034E-02, 8.4361E-02, 8.3272E-02,  &
        8.1899E-02, 8.0316E-02, 7.8614E-02, 7.6877E-02, 7.5151E-02,  &
        7.3461E-02, 7.1885E-02, 7.0510E-02, 6.9353E-02, 6.8352E-02,  &
        6.7484E-02, 6.6763E-02, 6.6218E-02, 6.5908E-02, 6.5843E-02,  &
        6.6050E-02, 6.6033E-02, 6.7196E-02, 6.9507E-02, 7.2385E-02,  &
        7.5103E-02, 7.7279E-02, 7.8414E-02, 7.8571E-02, 7.8530E-02,  &
        7.8530E-02 /)

real(f4), dimension(size(pbnd_101)) :: CH4_0 =  &  ! Default CH4 profile [ppm]
     (/ 2.8054E-02, 5.9711E-02, 1.2478E-01, 1.9952E-01, 2.3767E-01,  &
        2.5807E-01, 2.6910E-01, 2.7917E-01, 2.8537E-01, 2.9144E-01,  &
        3.0148E-01, 3.1228E-01, 3.2129E-01, 3.2810E-01, 3.3220E-01,  &
        3.3315E-01, 3.3026E-01, 3.2483E-01, 3.1810E-01, 3.1177E-01,  &
        3.0563E-01, 3.0137E-01, 3.0591E-01, 3.2537E-01, 3.5577E-01,  &
        4.1305E-01, 4.8083E-01, 5.6445E-01, 6.5112E-01, 7.4209E-01,  &
        8.2509E-01, 9.1022E-01, 9.8437E-01, 1.0619E+00, 1.1231E+00,  &
        1.1875E+00, 1.2369E+00, 1.2883E+00, 1.3257E+00, 1.3629E+00,  &
        1.3989E+00, 1.4357E+00, 1.4713E+00, 1.5066E+00, 1.5392E+00,  &
        1.5651E+00, 1.5922E+00, 1.6120E+00, 1.6321E+00, 1.6502E+00,  &
        1.6646E+00, 1.6796E+00, 1.6908E+00, 1.7002E+00, 1.7100E+00,  &
        1.7176E+00, 1.7248E+00, 1.7322E+00, 1.7389E+00, 1.7454E+00,  &
        1.7522E+00, 1.7583E+00, 1.7636E+00, 1.7692E+00, 1.7748E+00,  &
        1.7796E+00, 1.7846E+00, 1.7896E+00, 1.7943E+00, 1.7982E+00,  &
        1.8022E+00, 1.8063E+00, 1.8096E+00, 1.8121E+00, 1.8145E+00,  &
        1.8170E+00, 1.8192E+00, 1.8205E+00, 1.8218E+00, 1.8231E+00,  &
        1.8244E+00, 1.8255E+00, 1.8264E+00, 1.8274E+00, 1.8284E+00,  &
        1.8293E+00, 1.8299E+00, 1.8306E+00, 1.8312E+00, 1.8323E+00,  &
        1.8336E+00, 1.8350E+00, 1.8368E+00, 1.8391E+00, 1.8414E+00,  &
        1.8426E+00, 1.8431E+00, 1.8431E+00, 1.8431E+00, 1.8431E+00,  &
        1.8431E+00 /)


contains


subroutine init_PCRTM_invocation(PCRTM_static_input_dir, calc_rad_chan,  &
     enable_bT_chan, enable_Jacobian_chan, calc_bT_Jacobian, calc_tr_chan,  &
     molindx, scalfac, &
     PCRTM_STND, ATM_ABS_COEF, ICE_GRID, WAT_GRID, EOF_SOLUTION, ATM,  &
     RT_SOLUTION, K_M, K_PC, K_CH, TR_SOLUTION, nPCbnd, jacob_val)

use INIT_PCRTM

implicit none

character(len=120), intent(IN) :: PCRTM_static_input_dir  ! Length in PCRTM lib
logical(b4), intent(IN) :: calc_rad_chan, enable_bT_chan,  &
     enable_Jacobian_chan, calc_bT_Jacobian, calc_tr_chan, jacob_val
integer(i4), dimension(max_n_abs_molec), intent(IN) :: molindx
real(f4), dimension(max_n_abs_molec), intent(IN) :: scalfac
type(PCRTM_ATM_ABSORPTION_TYPE), intent(INOUT) :: atm_abs_coef
type(PCRTM_TR_SOLUTION_TYPE), allocatable, dimension(:), intent(INOUT) ::  &
     tr_solution
type(PCRTM_CLD_TABLE_DEF), intent(INOUT) :: ice_grid, wat_grid
type(PCRTM_EOF_SOLUTION_TYPE), allocatable, dimension(:), intent(INOUT) ::  &
     EOF_solution
type(PCRTM_RT_SOLUTION_TYPE), intent(INOUT) :: RT_solution
type(PCRTM_NM_JACOBIAN_TYPE), intent(INOUT) :: K_M
type(PCRTM_CH_JACOBIAN_TYPE), allocatable, dimension(:), intent(INOUT) :: K_CH
type(PCRTM_PC_JACOBIAN_TYPE), allocatable, dimension(:), intent(INOUT) :: K_PC
type(PCRTM_ATMOSPHERE_TYPE), intent(INOUT) :: atm
type(PCRTM_ATM_ABS_STRUCT_TYPE), intent(INOUT) :: PCRTM_STND
integer(i4), allocatable, dimension(:), intent(INOUT) :: nPCbnd

if (allocated(EOF_solution)) return  ! Nothing to appears to need initialized

! INSTANTIATE PCRTM LUTS AND INPUT OUTPUT SOLUTION VARIABLES
CALL PCRTM_INIT(PCRTM_static_input_dir, sensor_ID, JACOB_val,  &
     calc_rad_chan, enable_bT_chan, enable_Jacobian_chan,  &
     calc_bT_Jacobian, calc_tr_chan,  &
     PCRTM_STND, ATM_ABS_COEF, ICE_GRID, WAT_GRID,  &
     EOF_SOLUTION, ATM, RT_SOLUTION, K_M, K_PC, K_CH,  &
     TR_SOLUTION)
!                FRQCH)  ! FRQCH IS OPTIONAL INPUT

! ASSIGN ABSORBING MOLECULE INDEX AND SCALING FACTOR
PCRTM_STND%MOLINDX(1:max_n_abs_molec) = MOLINDX(1:max_n_abs_molec)
PCRTM_STND%SCALFAC(1:max_n_abs_molec) = SCALFAC(1:max_n_abs_molec)

! DEFINE THE NUMBER OF EOFS TO BE USED IN EACH BAND. THIS PART CAN BE NEGLECTED 
! IF ALL EOFS ARE GOING TO BE USED  
allocate(nPCbnd(SIZE(EOF_SOLUTION)))
nPCbnd(:) = 100

CALL TRUNCATE_PCRTM_EOF_JACOB(EOF_SOLUTION, JACOB_val, K_PC, nPCbnd)

end subroutine init_PCRTM_invocation


subroutine invoke_PCRTM(PCRTM_STND, ATM_ABS_COEF, ICE_GRID, WAT_GRID,  &
     EOF_SOLUTION, ATM, RT_SOLUTION, K_M, K_PC, K_CH,  &
     TR_SOLUTION, atm_1Dfields, cld_1Dfields, emis, p_sfc, satang, T_skin,  &
     newCO2, n_profile_Jacobians)

use CLEAR_PCRTM
use PCRTM_FORWARD_MODEL

implicit none

type(PCRTM_ATM_ABSORPTION_TYPE), intent(INOUT) :: atm_abs_coef
type(PCRTM_TR_SOLUTION_TYPE), allocatable, dimension(:), intent(INOUT) ::  &
     tr_solution
type(PCRTM_CLD_TABLE_DEF), intent(INOUT) :: ice_grid, wat_grid
type(PCRTM_EOF_SOLUTION_TYPE), allocatable, dimension(:), intent(INOUT) ::  &
     EOF_solution
type(PCRTM_RT_SOLUTION_TYPE), intent(INOUT) :: RT_solution
type(PCRTM_NM_JACOBIAN_TYPE), intent(INOUT) :: K_M
type(PCRTM_CH_JACOBIAN_TYPE), allocatable, dimension(:), intent(INOUT) :: K_CH
type(PCRTM_PC_JACOBIAN_TYPE), allocatable, dimension(:), intent(INOUT) :: K_PC
type(PCRTM_ATMOSPHERE_TYPE), intent(INOUT) :: atm
type(PCRTM_ATM_ABS_STRUCT_TYPE), intent(INOUT) :: PCRTM_STND
real(f8), dimension(:,:), intent(IN) :: atm_1Dfields, cld_1Dfields, emis
real(f4), intent(IN) :: T_skin, p_sfc, newCO2, satang
integer(i4), intent(IN) :: n_profile_Jacobians

! Prepare parameters/profiles for RTM run:

GEOMETRY%POBS = 0.005  ! [hPa] OBSERVATION PRESSURE
Geometry%psfc = p_sfc  ! [hPa]
Geometry%satang = satang  ! [deg]

RT_solution%emis(:) = emis(:,1)  ! For current sfctype only

Atm%Tskin = T_skin  ! [K]
Atm%Tlev(:) = atm_1Dfields(:,1)  ! [K]

Atm%VMR(:,1) = atm_1Dfields(:,2)  ! H2Ov [g/kg]
Atm%VMR(:,2) = CO2_0(:)*newCO2/385.85  ! CO2 [ppm]
Atm%VMR(:,3) = atm_1Dfields(:,3)  ! O3 [ppm]
Atm%VMR(:,4) = N2O_0(:)  ! N2O [ppm]
Atm%VMR(:,5) = CO_0(:)  ! CO [ppm]
Atm%VMR(:,6) = CH4_0(:)  ! CH4 [ppm]

Atm%pcld(:) = cld_1Dfields(:,1)
Atm%cld_flag(:) = int(cld_1Dfields(:,2))
Atm%taucld(:) = cld_1Dfields(:,3)
Atm%decld(:) = cld_1Dfields(:,4)

! Run RTM:
call PCRTM_FORWARD_RT(n_profile_Jacobians, GEOMETRY, ATM, ATM_ABS_COEF,  &
                      PCRTM_STND, ICE_GRID,  &
                      WAT_GRID, RT_SOLUTION, EOF_SOLUTION, K_M, K_PC, K_CH,  &
                      tr_SOLUTION)

end subroutine invoke_PCRTM


subroutine SFC_driver_for_PCRTM(n_i1fi, i1_flat_input, n_f8fi, f8_flat_input,  &
                                n_lev, n_atm_1Dfields, atm_1Dfields,  &
                                n_cld, n_cld_1Dfields, cld_1Dfields,  &
                                n_e1, n_e2, emis, PCRTM_static_input_dir,  &
                                n0_a1, n_a1, answer1, answer2)

implicit none

integer, intent(IN) :: n_i1fi, n_f8fi, n_lev, n_atm_1Dfields, n_cld,  &
     n_cld_1Dfields, n_e1, n_e2
integer(1), dimension(n_i1fi), intent(IN) :: i1_flat_input
real(f8), dimension(n_f8fi), intent(IN) :: f8_flat_input
real(f8), dimension(n_lev,n_atm_1Dfields), intent(IN) :: atm_1Dfields
real(f8), dimension(n_cld,n_cld_1Dfields), intent(IN) :: cld_1Dfields
real(f8), dimension(n_e1,n_e2), intent(IN) :: emis
character(len=120), intent(IN) :: PCRTM_static_input_dir  ! Length in PCRTM lib
integer, intent(IN) :: n0_a1
integer, intent(OUT) :: n_a1

real(f8), dimension(n0_a1), intent(INOUT) :: answer1, answer2

integer(i4), dimension(max_n_abs_molec) :: molindx
real(f4), dimension(max_n_abs_molec) :: scalfac
logical(b4) :: calc_tr_chan, calc_rad_chan, enable_Jacobian_chan,  &
     enable_bT_chan, calc_bT_Jacobian, calc_some_Jacobians, jacob_val
real(f4) :: T_skin, p_sfc, newCO2, satang
integer(i4) :: sfctype, n_profile_Jacobians
integer(i4), dimension(9) ::  &
     calc_these_Jacobians  ! H2O CO2 O3 N2O CO CH4 T Ts emis
integer(i4) :: im, ib, i1, ip

!=== Unpack some input fields:

ip = 1
molindx(1:max_n_abs_molec) = i1_flat_input(ip:ip+max_n_abs_molec-1)
ip = ip+max_n_abs_molec
calc_some_Jacobians = ( i1_flat_input(ip) == 1 ); ip = ip+1
calc_these_Jacobians(1:9) = i1_flat_input(ip:ip+9-1); ip = ip+9
enable_Jacobian_chan = ( i1_flat_input(ip) == 1 ); ip = ip+1
calc_bT_Jacobian = ( i1_flat_input(ip) == 1 ); ip = ip+1
enable_bT_chan = ( i1_flat_input(ip) == 1 ); ip = ip+1
calc_tr_chan= ( i1_flat_input(ip) == 1 ); ip = ip+1
calc_rad_chan= ( i1_flat_input(ip) == 1 ); ip = ip+1
sfctype = i1_flat_input(ip)

ip = 1
scalfac(1:max_n_abs_molec) = f8_flat_input(ip:ip+max_n_abs_molec-1)
ip = ip+max_n_abs_molec
T_skin = f8_flat_input(ip); ip = ip+1
p_sfc = f8_flat_input(ip); ip = ip+1
newCO2 = f8_flat_input(ip); ip = ip+1
satang = f8_flat_input(ip); ip = ip+1

  ! Number of profile Jacobians to compute; always zero for this application
n_profile_Jacobians = 0_i4

if (calc_some_Jacobians) then  ! Calculate one or more Jacobians

   jacob_val = .true.
   call init_PCRTM_invocation(PCRTM_static_input_dir, calc_rad_chan,  &
        enable_bT_chan, enable_Jacobian_chan, calc_bT_Jacobian, calc_tr_chan,  &
        molindx, scalfac, PCRTM_STND_Ja, ATM_ABS_COEF_Ja, ICE_GRID_Ja,  &
        WAT_GRID_Ja, EOF_SOLUTION_Ja, ATM_Ja, RT_SOLUTION_Ja, K_M_Ja,  &
        K_PC_Ja, K_CH_Ja, TR_SOLUTION_Ja, nPCbnd_Ja, jacob_val)

   call invoke_PCRTM(PCRTM_STND_Ja, ATM_ABS_COEF_Ja, ICE_GRID_Ja,  &
        WAT_GRID_Ja, EOF_SOLUTION_Ja, ATM_Ja, RT_SOLUTION_Ja, K_M_Ja,  &
        K_PC_Ja, K_CH_Ja, TR_SOLUTION_Ja, atm_1Dfields, cld_1Dfields, emis,  &
        p_sfc, satang, T_skin, newCO2, n_profile_Jacobians)

   if (calc_rad_chan) then
      !!!!@@@@ Done for JACOBIAN_CASE

      if (calc_these_Jacobians(8) == 1)  then
         ip = 0
         do ib=1,size(EOF_solution_Ja)
            do im=1,EOF_solution_Ja(ib)%nchbnd
               ip = ip+1
               answer1(ip) = EOF_solution_Ja(ib)%frqch(im)
               answer2(ip) = K_ch_Ja(ib)%R_ts(im)
            end do
         end do
         n_a1 = ip
!         print *, 'Jacobian_case (sfc_T)', size(EOF_solution_Ja), n_a1
      end if
                   
      if (calc_these_Jacobians(9) == 1)  then
         ip = 0
         do ib=1,size(EOF_solution_Ja)
            do im=1,EOF_solution_Ja(ib)%nchbnd
               ip = ip+1
               answer1(ip) = EOF_solution_Ja(ib)%frqch(im)
               answer2(ip) = K_ch_Ja(ib)%R_em(im)
            end do
         end do
         n_a1 = ip
!         print *, 'Jacobian_case (sfc_emis)', size(EOF_solution_Ja), n_a1
      end if

   end if

else   ! Forward model

   jacob_val = .false.
   call init_PCRTM_invocation(PCRTM_static_input_dir, calc_rad_chan,  &
        enable_bT_chan, enable_Jacobian_chan, calc_bT_Jacobian, calc_tr_chan,  &
        molindx, scalfac, PCRTM_STND, ATM_ABS_COEF, ICE_GRID, WAT_GRID,  &
        EOF_SOLUTION, ATM, RT_SOLUTION, K_M, K_PC, K_CH, TR_SOLUTION, nPCbnd,  &
        jacob_val)

   call invoke_PCRTM(PCRTM_STND, ATM_ABS_COEF, ICE_GRID, WAT_GRID,  &
        EOF_SOLUTION, ATM, RT_SOLUTION, K_M, K_PC, K_CH,  &
        TR_SOLUTION, atm_1Dfields, cld_1Dfields, emis, p_sfc, satang, T_skin,  &
        newCO2, n_profile_Jacobians)

   if (calc_rad_chan) then
      !!!!@@@@ Done for FORWARD_CASE
      ip = 0
      do ib=1,size(EOF_solution)
         do im=1,EOF_solution(ib)%nchbnd
            ip = ip+1
            answer1(ip) = EOF_solution(ib)%frqch(im)
            answer2(ip) = EOF_solution(ib)%radch(im)
         end do
      end do
      n_a1 = ip
!      print *, 'forward_case', size(EOF_solution), n_a1
   end if

end if
               
!@CALL CLEAR_PCRTM_LP( PCRTM_STND,    &
!@                       ICE_GRID,      &
!@                       WAT_GRID,      &
!@                       ATM_ABS_COEF,  &
!@                       ATM,           & 
!@                       RT_SOLUTION,   &
!@                       EOF_SOLUTION,  &
!@                       TR_SOLUTION )

!@IF (JACOB) THEN
!@     CALL CLEAR_PCRTM_JACOB( K_M,        &
!@                             K_PC,       &
!@                             K_CH )
!@END IF

end subroutine SFC_driver_for_PCRTM

end module SFC_driver_for_PCRTM_mod
