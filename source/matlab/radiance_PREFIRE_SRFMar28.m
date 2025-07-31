function [wv_center, rad_PREFIRE] = radiance_PREFIRE_SRFMar28(wn, radiance, ...
                                                              SRF, footprint_id)
% This code is to convert radiance at PCRTM wavenumber to radiance at PREFIRE channel wavelength 
%  
% input %%%%%%%
% radiance    unit is mWm-2/sr/cm-1
% SRF   spectral response info for PREFIRE
% updated on Mar.28, 2021 which has 3 dimensions of spectrum, channel, footprint
% footprint_id   the id_th footprint of SRF
%%%%%%%%%%%%%%%

%%%%%%% output
% wv_center  wavelength at the center of each channel, micron
% rad_PREFIRE radiance, unit is Wm-2/um/sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channel_wavelens(:,1) = SRF.channel_wavelen1_T(:,footprint_id);
channel_wavelens(:,2) = SRF.channel_wavelen2_T(:,footprint_id);

tmp = SRF.srf(footprint_id,:,:);
SRF_values = reshape(tmp,size(tmp,2), size(tmp,3))';
clear tmp;

wavelen = SRF.wavelen;

wv_center = 0.5*(channel_wavelens(:,1) + channel_wavelens(:,2));

nch = size(SRF_values, 2); 
dwv =wavelen(2) - wavelen(1);  % wavelength interval
 
nw = length(wn);
 
rad2 = radiance .* 1e-3 .* wn.^2/10000;  % convert from mWm-2/sr/cm-1 to Wm-2/sr/um
        
rad3 = interp1(10000./wn, rad2, wavelen);        
for ich =1:nch
        idw = find(wavelen>=channel_wavelens(ich,1) & wavelen<=channel_wavelens(ich,2));
        if sum(SRF_values(idw,ich))>0 & ~isnan(sum(rad3(idw)))
               
               rad_PREFIRE(ich) = nansum(rad3 .* SRF_values(:,ich) .* dwv)/ ...
                                  nansum(SRF_values(:,ich) .* dwv);
               
            else
                rad_PREFIRE(ich) = NaN;
            end
 end
