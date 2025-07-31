function [P, wv_center, Jacobian_PREFIRE] = Jacobian_PREFIRE_for_T_q_SRFMar28(wn, Jac, srf_name,footprint_id)
% This code is to convert Jacobian at PCRTM wavenumber to Jacobian at PREFIRE channel wavelength
% for T or q (h2o) at pressure level 

% input %%%%%%%
% Jac   Jacobian at pressure level from PCRTM simulation 
% srf_name     the file name of spectral response function for PREFIRE
% updated on Mar.28, 2021 which has 3 dimensions of spectrum, channel, footprint
% footprint_id   the id_th footprint of SRF
%%%%%%%%%%%%%%%

%%%%%%% output
% P      Pressure (hPa)
% wv_center  wavelength at the center of each channel, micron
% Jacobian_PREFIRE Jacobian for T or q at pressure level, unit is
%                  Wm-2/um/sr/K for T, and Wm-2/um/sr/(g/kg) for q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channel_wavelens(:,1)= ncread(srf_name,'channel_wavelen1',[footprint_id 1], [1 Inf])';

channel_wavelens(:,2) = ncread(srf_name,'channel_wavelen2',[footprint_id 1], [1 Inf])';

tmp = ncread(srf_name,'srf',[footprint_id 1 1],[1 Inf Inf]);
SRF = reshape(tmp,size(tmp,2), size(tmp,3))';
clear tmp;

wavelen = ncread(srf_name,'wavelen');

wv_center = 0.5*(channel_wavelens(:,1) + channel_wavelens(:,2));

nch = size(SRF,2); 
dwv =wavelen(2) - wavelen(1);  % wavelength interval
 
nw = length(wn);

 % PCRTM model level pressure in mb
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

 nlev = size(Jac,1);
 P = Pres(1:nlev);  
 
 Jacobian_PCRTM = Jac(:,1:nw)' .* 1e-3 .* wn.^2/10000;  % convert from mWm-2/sr/cm-1/X to Wm-2/sr/um/X (where x is K for T, and g/kg for h2o)
                           
 for ich =1:nch
        idw = find(wavelen>=channel_wavelens(ich,1) & wavelen<=channel_wavelens(ich,2));
        if nansum(SRF(idw,ich))>0
                for ilev=1:nlev
                    Jacobian_PCRTM2 = interp1(10000./wn, Jacobian_PCRTM(:,ilev), wavelen);
                    Jacobian_PREFIRE(ich,ilev) = nansum(Jacobian_PCRTM2 .* SRF(:, ich) .* dwv)/nansum(SRF(:, ich) .* dwv);
                end
        else
                Jacobian_PREFIRE(ich,1:nlev) = NaN;
        end
            
 end
               
      
