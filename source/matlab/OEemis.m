function [wv_emis,wv_rad,wn_emis,wn_rad,x,y,S_aposteriori,x_op,y_op,S_op,...
    K_op,dgf,converged,iter,chi2,pvalue,c2,ec] = ...
    OEemis(pcrtmpath, ancpath,footprint_id,valid_rad,valid_emis,...
    rad_obs,lat,lon,sfctype,pres,Ps,q,ozone,T,Ts, prior_emis_profile, SRF, ...
    MODTRAN_default_profile, emis_id2, ec)
% this function prepares the initial guess, a priori mean, a priori
% covariance, and error covariance for the surface emissivity retrieval test
%% Input
% pcrtmpath: filepath of the PCRTM package
% ancpath: ancillary data folder
% footprint_id: from 1 to 8
% valid_rad: PREFIRE channels for radiances, size 14*1
% valid_emis: PREFIRE channels for surface emissivity, size 14*1
% rad_obs: observed radiance spectrum, size 63 * 1, unit W m^-2 sr^-1 um^-1
% lat: latitude, from -90 to 90
% lon: longitude, from  0 to 360
% sfctype: surface type in ./PCRTM_V3.4/data_made_by_Chen/EMISSIVITY_id2
% pres: pressure levels, unit hPa
% Ps: surface pressure, unit hPa
% q: Specific humidity profile, unit: g/kg
% ozone: Ozone profile, unit: g/kg
% T: Temperature profile, unit: K
% Ts: Surface temperature, unit: K
% prior_emis_profile: prior emissivity profile information: -
% SRF: structure array containing various SRF-related fields
% MODTRAN_default_profile: structure array containing std MODTRAN atmos profiles
% emis_id2: PCRTM wavenumbers and stock emissivities for sensor id2

%% Output
% wv_emis: central wavelength of PREFIRE channels for surface emissivity, size 14*1, unit um
% wv_rad: central wavelength of PREFIRE channels for radiance, size 14*1, unit um
% wn_emis: central wavenumber of PREFIRE channels for surface emissivity, size 14*1, unit cm^-1
% wn_rad: central wavenumber of PREFIRE channels for radiance, size 14*1, unit cm^-1
% x: surface emissivity retrievals at each iteration step
% y: PCRTM-generated radiance based on surface emissivity retrieval at each iteration step
% S_aposteriori: A posteriori covariance matrix at each iteration step
% x_op: optimal estimate of surface emissivity
% y_op: PCRTM-generated radiance based on optimal estimated surface emissivity
% S_op: A posteriori covarinace matrix of optimal estimated surface emissivity
% K_op: Jacobian of optimal estimated surface emissivity
% dgf: degree of freedom for the optimal estimated surface emissivity
% converged: false / true, whether the retrieval is converged
% iter: iterations when the retrieval is converged
% chi2: Chisquare values, cost function in x-space and y-space
% pvalue: corresponding p-values of the chi-square test (chi2)
% c2: the nonlinearity parameter, Section 5.1 (Rodger 2000)

%% test of Optimal Estimation Method
maxIter = 100;                                                            % maximum number of iterations
maxTime = 1e7;                                                            % maximum runtime
gammaFactor = ones(1,maxIter);                                            % gamma factor
gammaFactor(1:6) = [1000,300,100,30,10,3];
convergenceFactor = 10;

nemis = length(valid_emis);
nrad = length(valid_rad);

wn = (SRF.channel_wavenum1_T + SRF.channel_wavenum2_T)./2;
wn_emis = wn(valid_emis);
wn_rad = wn(valid_rad);
wv = (SRF.channel_wavelen1_T + SRF.channel_wavelen2_T)./2;  % central wavelength
wv_emis = wv(valid_emis);
wv_rad = wv(valid_rad);

%% prior mean covariance
xa = 0.95 * ones(nemis,1);                                                % a priori state value equal to the first guess                   

%---------Determine prior mean and variance for surface emissivity---------%
emis_apriori = 4 * cov(prior_emis_profile(:,valid_emis)); 

% % diagonal covariance matrix: 0.0225 on its diagonal
% Sa(1:nemis,1:nemis) = diag(0.1.^2 * ones(nemis,1));

% A priori covariance of surface emissivity has off-diagonal variables
Sa = 0.5 * emis_apriori;         % decrease off-diagonal correlation coefficients 
for ich = 1:nemis
   Sa(ich,ich) = emis_apriori(ich,ich);      % diagonal values: correlation coefficient 1    
end

%% Measurement error covariance
% Noise equivalent spectral radiance for PREFIRE
Se = diag(SRF.NEDR_T(valid_rad).^2);  
Sy = Se;                      % error contains only measurement noise (no model error)

%% Observation series
y_obs = rad_obs(valid_rad);                                                
y_obs = reshape(y_obs,nrad,1);

%% Initial Guess of the state vector: emissivity
x0 = xa;
    
%% Call the retrieval function
[x, y, S_aposteriori, x_op, y_op, S_op, K_op, dgf, converged, iter, chi2, ...
 pvalue, c2, ec] = retrievalPCRTM(pcrtmpath, ancpath, x0, xa, Sa, y_obs, Sy, ...
                gammaFactor, convergenceFactor, maxIter, maxTime, valid_rad, ...
                valid_emis, lat, lon, sfctype, pres, Ps, q, ozone, T, Ts, ...
                footprint_id, SRF, MODTRAN_default_profile, emis_id2, ec);
end
