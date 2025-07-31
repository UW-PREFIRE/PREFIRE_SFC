function emis_new = expandEmis(emis0, SRF, footprint_id, valid_emis)
% This code is based on changeEmis.m to expand the emissivity at 14
% retrived channels to all 63 PREFIRE channels

% emis0: surface emissivity at 14 selected PREFIRE channels
% wn0: wavenumber of the corresponding PREFIRE spectrum
%% 

% Updated SRFs on March 28, 2021 
wn_min = SRF.channel_wavenum1_T(:,footprint_id);
wn_max = SRF.channel_wavenum2_T(:,footprint_id);

wn_tot = cat(2,wn_min,wn_max);

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

emis_new(1:temp(1)-1) = emis(1);
emis_new(temp) =  emis;
emis_new(temp(end)+1:length(wn_min)) = emis(end);

