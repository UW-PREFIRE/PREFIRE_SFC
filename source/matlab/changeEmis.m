function [emis_new] = changeEmis(emis_id2, sfctype, emis0, SRF, ...
                                 footprint_id, valid_emis)
% The PCRTM simulator changes surface emissivity through surface type
% This function is developed to change corresponding emissivity column in
% "EMISSIVITY_id2"

% sfctype: the surface type of one grid box, integer from 1 to 18
% emis0: surface emissivity at 14 selected PREFIRE channels
% wn0: wavenumber of the corresponding PREFIRE spectrum
%% 

wn_tot = cat(2, SRF.channel_wavenum1_T(:,footprint_id), ...
             SRF.channel_wavenum2_T(:,footprint_id));

temp = valid_emis(1) : valid_emis(end);
wn_temp = wn_tot(temp,:);

% select out the index of non-valid PREFIRE channels within the limits
id_nonvalid = temp(~ismember(temp, valid_emis));

% fill in nonvalid PREFIRE channels within the limits 
emis = NaN(length(temp),1);
emis(valid_emis-valid_emis(1)+1) = emis0;

num = 1;
while num <= length(id_nonvalid)
   a = id_nonvalid(num);
   b = a - valid_emis(1)+1;
   ii = 1;
   while ismember(a+ii,id_nonvalid)
      ii = ii + 1; 
   end
   
   if mod(ii,2)==1
       if ii == 1
           emis(b) = (emis(b-1) + emis(b+1))./2;          
       else
           for jj = 1:((ii-1)/2)
               emis(b+jj-1) = emis(b+jj-2);
               emis(b+ii-jj) = emis(b+ii-jj+1);
           end
           emis(b+jj) = (emis(b+jj-1) + emis(b+jj+1))./2;
       end
       
   elseif mod(ii,2)==0
       for jj = 1:(ii/2)
           emis(b+jj-1) = emis(b+jj-2);
           emis(b+ii-jj) = emis(b+ii-jj+1);
       end      
   end
   
   num = num + ii;
end

%% Prepare the new emissivity of size 740*1
emis_new = zeros(length(emis_id2.wn_new),1);

% fill in the wavenumbers within the limits
for i = 1:length(temp)
    idw = find(emis_id2.wn_new >= wn_temp(i,1) & emis_id2.wn_new <= wn_temp(i,2));
    emis_new(idw,1) = emis(i);
end

% fill in non-valid wavenumbers out of the limits
id_less = find(emis_id2.wn_new < min(wn_temp(:,1)));
id_larger = find(emis_id2.wn_new > max(wn_temp(:,2)));
emis_new(id_less,1) = emis_new(id_less(end)+1,1);
emis_new(id_larger,1) = emis_new(id_larger(1)-1,1);
