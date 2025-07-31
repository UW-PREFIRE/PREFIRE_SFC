function radiance = forwardPCRTM(pcrtmpath,ancpath,valid_rad,valid_emis,x, ...
                                 lat,lon,sfctype,pres,Ps,q,ozone,TT,Ts, ...
                                 footprint_id, SRF, MODTRAN_default_profile, ...
                                 emis_id2)

%%%%%%%%%%%%%%%%%%%%% Configuration %%%%%%%%%%%%%%%%%%%%
% set parameters
VZA=0;
ncol=50;
overlap = 3; 
co2 = 380;  
ATMno = 6;   
Plen = 5421;  % number of PCRTM channels for sensor_id=2
max_n_abs_molec = 6;  % maximum number of absorbing molecules
%%%%%%%%%%%%%%%%%%% end of configuration%%%%%%%%%%%%%%%   

lev = double(pres); nlev = length(lev);
                
% Variable x : Emissivity
nemis = length(valid_emis);
emis0 = x(1:nemis);

nboxes =1;                                      
%@emis(nboxes) = 1.0; % a default value, not used in this case

o3 = ozone;
cc = zeros(1,nlev);
clwc = zeros(1,nlev);
ciwc = zeros(1,nlev);
    
% print out central latitude & longitude to make sure the region is right
lato(nboxes) = lat;
lono(nboxes) = lon;  

profiles_for_RTM = make_profiles(nboxes,lato,lono, nlev, ncol, overlap, clwc, ciwc, cc, lev, TT, q, o3, Ps, Ts, ATMno, sfctype, co2, pcrtmpath, VZA, MODTRAN_default_profile);

%%%%%%%%%%%%%%%%%%%%%% Change the EMISSIVITY_id2 %%%%%%%%%%%%%%%%%%%%
emis_new = changeEmis(emis_id2, sfctype, emis0, SRF, footprint_id, ...
                      valid_emis);

%% run PCRTM to calculate radiance

molindx = [2 2 2 2 2 2];  % H2O CO2 O3 N2O CO CH4
scalfac = [1.0, 1.0135, 1.0, 1.027, 1.0, 1.145];  % H2O CO2 O3 N2O CO CH4
calc_some_Jacobians = false;
calc_these_Jacobians = [0 0 0 0 0 0 0 0 0];  % H2O CO2 O3 N2O CO CH4 T Ts emis
enable_Jacobian_chan = false;
calc_bT_Jacobian = false;
enable_bT_chan = false;
calc_tr_chan = false;
calc_rad_chan = true;
satang = 0.;

% Preallocate some arrays:
tmp_i1_flat = zeros(max_n_abs_molec+9+7, 1);
tmp_f8_flat = zeros(max_n_abs_molec+4, 1);
atm_1Dfields = zeros(101,3);
max_n_clds = 101; cld_1Dfields = zeros(max_n_clds,5);

for ibox=1:nboxes

   for icol=1:length(profiles_for_RTM(ibox).col)

      % Construct input array of 8-bit (1-byte) integers:
      ip = 1;
      tmp_i1_flat(ip:ip+max_n_abs_molec-1) = molindx(1:max_n_abs_molec);
      ip = ip+max_n_abs_molec;
      tmp_i1_flat(ip) = calc_some_Jacobians; ip = ip+1;
      tmp_i1_flat(ip:ip+9-1) = calc_these_Jacobians(1:9); ip = ip+9;
      tmp_i1_flat(ip) = enable_Jacobian_chan; ip = ip+1;
      tmp_i1_flat(ip) = calc_bT_Jacobian; ip = ip+1;
      tmp_i1_flat(ip) = enable_bT_chan; ip = ip+1;
      tmp_i1_flat(ip) = calc_tr_chan; ip = ip+1;
      tmp_i1_flat(ip) = calc_rad_chan; ip = ip+1;
      tmp_i1_flat(ip) = profiles_for_RTM(ibox).col(icol).sfctype;
      i1_flat_input = int8(tmp_i1_flat);

      % Construct input array of 8-byte floating point values:
      ip = 1;
      tmp_f8_flat(ip:ip+max_n_abs_molec-1) = scalfac(1:max_n_abs_molec);
      ip = ip+max_n_abs_molec;
      tmp_f8_flat(ip) = profiles_for_RTM(ibox).col(icol).Ts; ip = ip+1;
      tmp_f8_flat(ip) = profiles_for_RTM(ibox).col(icol).ps; ip = ip+1;
      tmp_f8_flat(ip) = profiles_for_RTM(ibox).col(icol).CO2_v; ip = ip+1;
      tmp_f8_flat(ip) = satang; ip = ip+1;

      atm_1Dfields(:,1) = profiles_for_RTM(ibox).col(icol).T_prf;  % air T [K]
      atm_1Dfields(:,2) = profiles_for_RTM(ibox).col(icol).H2Ov_prf;  % H2Ov [g/kg]
      atm_1Dfields(:,3) = profiles_for_RTM(ibox).col(icol).O3_prf;  % O3 [ppm]

      if profiles_for_RTM(ibox).col(icol).has_clouds
         cldnum = profiles_for_RTM(ibox).col(icol).cldnum;
         cld_1Dfields(1:cldnum,1) = profiles_for_RTM(ibox).col(icol).cldpres(:);
         cld_1Dfields(1:cldnum,2) = profiles_for_RTM(ibox).col(icol).cldopt(:);
         cld_1Dfields(1:cldnum,3) = profiles_for_RTM(ibox).col(icol).cldde(:);
         cld_1Dfields(1:cldnum,4) = profiles_for_RTM(ibox).col(icol).cld_id(:);
         cld_1Dfields(1:cldnum,5) = profiles_for_RTM(ibox).col(icol).cldphase(:);
      end

      PCRTM_static_input_dir = fullfile(pcrtmpath, 'INPUTDIR');
      [wn, tmprad] = PCRTM_MEX_interface(i1_flat_input, tmp_f8_flat, ...
                                   atm_1Dfields, cld_1Dfields, emis_new, ...
                                   PCRTM_static_input_dir);

      if icol ~= 1
         disp('Unknown how to handle this icol')
         exit(1)
      end

      % clear sky spectra
      [wv_center, rad_out] = radiance_PREFIRE_SRFMar28(wn, tmprad, ...
                                                       SRF, footprint_id);
      radiance(:,ibox) = rad_out;                      % clear-sky radiance in Wm-2/um/sr

   end
end

% radiance of valid PREFIRE channels
radiance = radiance(valid_rad,:);  
