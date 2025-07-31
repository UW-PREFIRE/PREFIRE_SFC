function [wv_center, Jacobian_PREFIRE] = ...
              Jacobian_PREFIRE_for_Ts_emis_SRFMar28(wn, Jac, SRF, footprint_id)
% This code is to convert Jacobian at PCRTM wavenumber to Jacobian at PREFIRE channel wavelength 
% for Ts or surface emissivity at surface level 

% input %%%%%%%
% Jac    Jacobian at surface level from PCRTM simulation
% SRF   spectral response info for PREFIRE
% updated on Mar.28, 2021 which has 3 dimensions of spectrum, channel, footprint
% footprint_id   the id_th footprint of SRF
%%%%%%%%%%%%%%%

%%%%%%% output
% wv_center  wavelength at the center of each channel, micron
% Jacobian_PREFIRE Jacobian for Ts or surface emissivity at surface level, unit is
%                  Wm-2/um/sr/K for Ts, and Wm-2/um/sr/1 for surface
%                  emissivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channel_wavelens(:,1) = SRF.channel_wavelen1_T(:,footprint_id);
channel_wavelens(:,2) = SRF.channel_wavelen2_T(:,footprint_id);

tmp = SRF.srf(footprint_id,:,:);
SRF_values = reshape(tmp,size(tmp,2), size(tmp,3))';
clear tmp;

wavelen = SRF.wavelen;

wv_center = 0.5*(channel_wavelens(:,1) + channel_wavelens(:,2));

nch = size(SRF_values, 2); 
dwv = wavelen(2) - wavelen(1);  % wavelength interval

nw = length(wn);
  
% convert from mWm-2/sr/cm-1/X to Wm-2/sr/um/X (where x is K for Ts, and 1 for surface emissivity)
Jacobian = Jac .* 1e-3 .* wn.^2/10000;  
Jacobian2 = interp1(10000./wn, Jacobian, wavelen);

for ich =1:nch
   idw = find(wavelen>=channel_wavelens(ich,1) & wavelen<=channel_wavelens(ich,2));
   if sum(SRF_values(idw,ich))>0 %& ~isnan(sum(Jacobian2(idw)))
      Jacobian_PREFIRE(ich) = nansum(Jacobian2 .* SRF_values(:,ich) .* dwv)/ ...
                              nansum(SRF_values(:,ich) .* dwv);               
        else
               Jacobian_PREFIRE(ich) = NaN;
        end
 end
