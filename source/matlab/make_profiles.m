function [profiles_for_RTM] = make_profiles( ...
      nboxes, lato,lono, nlev, ncol, overlap, WCT, ICT, cc, P, TT, q, o3, Ps, Ts, ATMno, sfctype, newco2, pcrtmpath, VZA, MODTRAN_default_profile)

% nboxes   the number of boxes
% lato     latitude of each boxes
% lono     longitude of each boxes
% nlev     the number of layers, all boxes has same levels
% ncol     the number of subcolumn  
% overlap =1 maximum overlap
% overlap =2 random overlap
% overlap =3 maximum-random overlap
% WCT      water content in kg/kg from top to bottom
% ICT      ice content in kg/kg from top to bottom
% cc       cloud fraction profile, 0-1, from top to bottom
% P        pressure profile in hPa from top to bottom
% TT       temperature profile in k, from top to bottom
% q        specific humidity in g/kg, from top to bottom
% o3       ozone mass mixing ratio in g/kg, from top to bottom
% Ps       surface pressure in hPa
% Ts       surface temperature in k
% ATMno    default model for the upper atmosphere
% sfctype  the surface type of each boxes
% emis        the surface emissivity
% if sfctype<0, use emissivity defined by emis
% VZA      viewing zentih angle
% MODTRAN_default_profile: structure array containing std MODTRAN atmos profiles

% PCRTM fixed level pressure in mb (total 101 levels)
% T,q, o3 profiles need to be interpolated into these fixed levels
% cloud profiles not need to be interpolated, but need level_id which shows
% cloud are between Pres(level_id) and Pres(level_id +1)
Pres = [    0.0050,    0.0161,    0.0384,    0.0769,    0.1370,    ...
            0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    ...
            1.2972,    1.6872,    2.1526,    2.7009,    3.3398,    ...
            4.0770,    4.9204,    5.8776,    6.9567,    8.1655,    ...
            9.5119,   11.0038,   12.6492,   14.4559,   16.4318,    ...
            18.5847,   20.9224,   23.4526,   26.1829,   29.1210,   ...
            32.2744,   35.6505,   39.2566,   43.1001,   47.1882,   ...
            51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   ...
            77.2396,   83.2310,   89.5204,   96.1138,  103.0172,   ...
            110.2366,  117.7775,  125.6456,  133.8462,  142.3848,  ...
            151.2664,  160.4959,  170.0784,  180.0183,  190.3203,  ...
            200.9887,  212.0277,  223.4415,  235.2338,  247.4085,  ...
            259.9691,  272.9191,  286.2617,  300.0000,  314.1369,  ...
            328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  ...
            407.4738,  424.4698,  441.8819,  459.7118,  477.9607,  ...
            496.6298,  515.7200,  535.2322,  555.1669,  575.5248,  ...
            596.3062,  617.5112,  639.1398,  661.1920,  683.6673,  ...
            706.5654,  729.8857,  753.6275,  777.7897,  802.3714,  ...
            827.3713,  852.7880,  878.6201,  904.8659,  931.5236,  ...
            958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170,  ...
            1100.0000 ]; 

R = 287.05;
g = 9.806;

% coefficents from S-C Ou, K-N. Liou, Atmospheric Research
% 35(1995):127-138.
% for computing ice cloud effective size
c0 = 326.3;
c1 = 12.42;
c2 = 0.197;
c3 = 0.0012;

%  cloud assignment to each sub column wthin a gridbox
[unique_col_frac ucol_num ucol_num_same] = get_subcolumn_frac_v2(nboxes, nlev, ncol, cc, overlap);  

% convert ozone from g/kg to ppmv
o3 = o3 *1e3/47.998*28.96; 

newT = zeros(1,101);
newh2o = zeros(1,101);
newozone = zeros(1,101);

profiles_for_RTM = [struct];

% use default profiles for Pres smaller than the smallest pressure of each
% box, as the value is identical for each box, we use P(1) as the smallest pressure
indU = find (P(1) >= Pres);

% !!!Changed by Yan
if isempty(indU)
    ii = 0;
else
    ii = indU(end);
end
% !!!

def_P = MODTRAN_default_profile.Pres(:, ATMno);
def_T = MODTRAN_default_profile.T_z(:, ATMno +1);
def_h2o = MODTRAN_default_profile.h2o(:, ATMno);
def_o3 = MODTRAN_default_profile.o3(:, ATMno);
% convert default h2o from ppmv to g/kg
def_h2o  = def_h2o *1e-3 *18/28.96;


% parameters for default models in upper atmosphere
% for ilev = 1: indU(end)
for ilev = 1:1:ii
        newT(ilev) = interp1(log(def_P), def_T, log(Pres(ilev)));
        newh2o(ilev) = interp1(def_P, def_h2o, Pres(ilev));
        newozone(ilev) = interp1(def_P, def_o3, Pres(ilev));
end


endsign=0;

for ibox =1:nboxes

    profiles_for_RTM(ibox).col = [struct];
    profiles_for_RTM(ibox).has_clouds = false;  % Initial assumption
    profiles_for_RTM(ibox).clearsky_icol = 0;  % Initial assumption

    %  set a minimum specific humidity as 1.0e-6 to avoid spcific humidity is zero    
    idw = find(q(ibox,:)<=0.0);
    q(ibox,idw) = 1.0e-6;
    % get T, q, o3 profiles for lower atmosphere
%     for ilev =indU(end) + 1:length(Pres)
    for ilev = (ii+1):length(Pres)
           % for Pres below the bottom pressure, T, q, and o3 are set to the value at bottom pressure        
         if (Pres(ilev)>= P(end))
             Pp = P(end);
         elseif (Pres(ilev)>= P(1) & Pres(ilev)< P(end))
             Pp = Pres(ilev);
         end
             
         newT(ilev) = interp1(log(P), TT(ibox,:), log(Pp));
         newh2o(ilev) = exp(interp1(log(P), log(q(ibox,:)), log(Pp)));         
         newozone(ilev) = exp(interp1(log(P), log(o3(ibox,:)), log(Pp)));  
    end
    
    % set a minimum cloud cover (0.001) to decide clear or cloudy over the box
    idc = find(cc(ibox,:)>0.001);

    % convert pressure level to altitude level
    Z = zeros(nlev,1);
         
    for ilev = length(P(:))-1:-1:1
	        Ttmp = (TT(ibox,ilev) + TT(ibox,ilev+1))/2;
            scaleH = R*Ttmp/g/1000;
            dz = scaleH*log(P(ilev+1)/P(ilev));
            Z(ilev) = Z(ilev+1) + dz; 
    end
      
    sub_opt(1:nlev) = 0.0;
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    % First, set cloud parameter over each gridbox
    
    % set cloud pressure as the midlevel pressure of the highest level, from Klein S.A. and Jakob C., 1999    
    cldpres_layer = exp(0.5*(log(P(1:end-1)) + log(P(2:end)))); 
    cldpres_layer(nlev) = P(nlev);
    % replace the last cloud layer pressure with surface pressure
    % May 27, 2012
    idp = find(P>=Ps(ibox));
    if ~isempty(idp)
            cldpres_layer(idp(1)-1) = Ps(ibox);
    end
   
    % set a default cloud effective size and cloud type on each layer
    cldde_layer(1:nlev) = 40.0;
    cldphase_layer(1:nlev) = 1;
         
    % get the corresponding id of P(ibox,:) to Pres
    %  P(ibox,i)>Pres(l) & P(ibox,i)>Pres(l) 
    %  this is for assigning cloud to the nearest level of atmosphere ( the fixed 101 levels in PCRTM)
    clear idp;
    for i=1 : length(cldpres_layer)
             
             itmp = find(cldpres_layer(i) < Pres);
             if ~isempty(itmp)
              idp(i) = itmp(1);
             end
    end
         
    WP(ibox) =0.0; 
    for i =1: length(idc)
             
               ilev = idc(i);
               %############################################################
               %compute cirrus cloud effective size from cloud temperature,
               % cite from S-C Ou, K-N. Liou, Atmospheric Research
               % 35(1995):127-138.
               % temperature of ice clouds is the atmosphere temperature
               % cirrus cloud effective size are fitted in the range from - 20 to - 60C 
               % here the range are reduced for most cirrsus temperature
               % is in this range
               tcld = TT(ibox, ilev) - 273.16;  % temperature of ice clouds
               if tcld<-50
                          tcld = -50;  % set a minimum 
               end
               if tcld>-25
                         tcld = -25;  % set a maximum
               end
                      
               cldde_ice(ilev) = c0 + c1 * tcld + c2 * tcld^2 + c3 * tcld^3 ;
                               
               
               cldde_liq(ilev) =20; % cloud effective size for water clouds, diameter in PCRTM
               
               % cloud phase and cloud effective size  in  micron
               if cldpres_layer(ilev) <=440
                         
                         cldphase_layer(ilev) = 1;  % cirrus cloud
                       
                         cldde_layer(ilev) = cldde_ice(ilev);
                     
                      
               else
                         cldphase_layer(ilev) = 2;  % water cloud
                         cldde_layer(ilev) = cldde_liq(ilev);
               end
               
              
               if ICT(ibox, ilev)>1e-10
                        % compute ice cloud optical depth from Ebert and Curry (1992, J. Geophys. Res., 
                        %    vol. 97, pp. 3831-3836.
                        % the same as GFDL
                        qi = ICT(ibox, ilev)/ TT(ibox,ilev) *P(ilev)*100/R *1e3;  %change ice water content from kg/kg to g/m^3  
                        ice_opt =  (0.003448 + 2.431*2/cldde_ice(ilev))*qi/cc(ibox,ilev)*(Z(ilev)-Z(min(ilev+1,nlev))) *1e3; % * 0.5 * (P(ilev)+P(max(ilev-1,1)))/g *1e4;
                    
               else
                        ice_opt =0;
                        qi = 0;
               end
               
             
               qw = WCT(ibox, ilev)/ TT(ibox,ilev) *P(ilev)*100/R *1e3;  %change liquid water content from kg/kg to g/m^3
               % cloud optical depth for water clouds, klein S. A. et al. (1999)                    
               % water_opt = 0.15893 * qw/cc(ibox,ilev) * 0.5 * (P(ilev)+P(max(ilev-1,1)))/g *1e4;
               
               % from ECMWF technical report                   
               water_opt = 3 * qw/cldde_liq(ilev)/cc(ibox,ilev)*(Z(ilev)-Z(min(ilev+1,nlev))) *1e3; % * 0.5 * (P(ilev)+P(max(ilev-1,1)))/g *1e4;
               sub_opt(ilev) = ice_opt + water_opt; 
               
                % cloud water path in g/m^2     
               WP(ibox) = WP(ibox) + (qi + qw) *(Z(ilev)-Z(min(ilev+1,nlev))) *1e3; %* 0.5 * (P(ilev)+P(max(ilev-1,1)))/g *1e4;
            
    end
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %  Now, set cloud parameters over each sub-column
    
    
    % prepare all sky profiles of each sub-column for PCRTM
    % only input one profile if some profiles are the same
    % ucol_num(ibox) is the number of unique profiles  
    for icol =1:ucol_num(ibox)  
        
            % find cloud levels over unique sub-column
            idx = find (unique_col_frac(ibox,icol,:) ==1);
            
            % If this sub-column has clouds
            if ~isempty(idx)
                
                cldnum = length(idx);  % number of cloud layers
                
                for ic =1: cldnum
                    
                     % cloud top pressure in hPa                   
                   
                     cldpres(ic) = cldpres_layer(idx(ic));
                     
                     cldphase(ic) = cldphase_layer(idx(ic));  
                     
                     cldde(ic) = cldde_layer(idx(ic));
                     cld_id(ic) = idp(idx(ic));
                     
                     % cloud optical depth
                     cldopt(ic) = sub_opt(idx(ic));       
                      
                end % cloud layers

               profiles_for_RTM(ibox).has_clouds = true;
               profiles_for_RTM(ibox).col(icol).has_clouds = true;
               profiles_for_RTM(ibox).col(icol).ps = Ps(ibox);
               profiles_for_RTM(ibox).col(icol).Ts = Ts(ibox);
               profiles_for_RTM(ibox).col(icol).emis = 1.;  % Dummy value
               profiles_for_RTM(ibox).col(icol).CO2_v = newco2;
               profiles_for_RTM(ibox).col(icol).sfctype = sfctype(ibox);
               profiles_for_RTM(ibox).col(icol).p_prf = Pres;
               profiles_for_RTM(ibox).col(icol).T_prf = newT;
               profiles_for_RTM(ibox).col(icol).H2Ov_prf = newh2o;
               profiles_for_RTM(ibox).col(icol).O3_prf = newozone;
               profiles_for_RTM(ibox).col(icol).cldpres = cldpres;
               profiles_for_RTM(ibox).col(icol).cldopt = cldopt;
               profiles_for_RTM(ibox).col(icol).cldde = cldde;
               profiles_for_RTM(ibox).col(icol).cldnum = cldnum;
               profiles_for_RTM(ibox).col(icol).cld_id = cld_id;
               profiles_for_RTM(ibox).col(icol).cldphase = cldphase;

            else
               profiles_for_RTM(ibox).clearsky_icol = icol;
               profiles_for_RTM(ibox).col(icol).has_clouds = false;
               profiles_for_RTM(ibox).col(icol).ps = Ps(ibox);
               profiles_for_RTM(ibox).col(icol).Ts = Ts(ibox);
               profiles_for_RTM(ibox).col(icol).emis = 1.;  % Dummy value
               profiles_for_RTM(ibox).col(icol).CO2_v = newco2;
               profiles_for_RTM(ibox).col(icol).sfctype = sfctype(ibox);
               profiles_for_RTM(ibox).col(icol).p_prf = Pres;
               profiles_for_RTM(ibox).col(icol).T_prf = newT;
               profiles_for_RTM(ibox).col(icol).H2Ov_prf= newh2o;
               profiles_for_RTM(ibox).col(icol).O3_prf = newozone;
            end
                
    end % unique cloud profile
         
   % prepare clear sky profile of each gridbox for PCRTM, if needed

    if profiles_for_RTM(ibox).clearsky_icol == 0
       icol = ucol_num(ibox)+1;
       profiles_for_RTM(ibox).clearsky_icol = icol;
       profiles_for_RTM(ibox).col(icol).has_clouds = false;
       profiles_for_RTM(ibox).col(icol).ps = Ps(ibox);
       profiles_for_RTM(ibox).col(icol).Ts = Ts(ibox);
       profiles_for_RTM(ibox).col(icol).emis = 1.;  % Dummy value
       profiles_for_RTM(ibox).col(icol).CO2_v = newco2;
       profiles_for_RTM(ibox).col(icol).sfctype = sfctype(ibox);
       profiles_for_RTM(ibox).col(icol).p_prf = Pres;
       profiles_for_RTM(ibox).col(icol).T_prf = newT;
       profiles_for_RTM(ibox).col(icol).H2Ov_prf= newh2o;
       profiles_for_RTM(ibox).col(icol).O3_prf = newozone;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
