function [x,y,S_aposteriori,x_op,y_op,S_op,K_op,dgf,converged,iter,chi2,pvalue,c2,ec] = ...
    retrievalPCRTM(pcrtmpath,ancpath,x0,xa,Sa,y_obs,Sy,gammaFactor,convergenceFactor,maxIter,~,...
    valid_rad,valid_emis,lat,lon,sfctype,pres,Ps,q,ozone,T,Ts,footprint_id, ...
    SRF, MODTRAN_default_profile, emis_id2, ec)
% This function performs optimal estimation retrieval
%% Input
% pcrtmpath: filepath of the PCRTM package
% ancpath: ancillary data folder
% x0: initial guess of surface emissivity
% xa: a priori mean of surface emissivity
% Sa: a priori covariance of surface emissivity
% y_obs: observed radiance on valid PREFIRE channels, size 14*1, unit W m^-2 sr^-1 um^-1
% Sy: measurement error covariance (no model error currently)
% gammaFactor: a tuning parameter vector in order to stablize the retrieval process
% convergenceFactor: a parameter used for convergence criterion 
% maxIter: maximum number of iterations
% '~' is preserved for maximum running time (now deleted)
% valid_rad: PREFIRE channels for radiances, size 14*1
% valid_emis: PREFIRE channels for surface emissivity, size 14*1
% lat: latitude, from -90 to 90
% lon: longitude, from  0 to 360
% sfctype: surface type in ./PCRTM_V3.4/data_made_by_Chen/EMISSIVITY_id2
% pres: pressure levels, unit hPa
% Ps: surface pressure, unit hPa
% q: Specific humidity profile, unit: g/kg
% ozone: Ozone profile, unit: g/kg
% T: Temperature profile, unit: K
% Ts: Surface temperature, unit: K
% footprint_id: from 1 to 8
% SRF: structure array containing various SRF-related fields
% MODTRAN_default_profile: structure array containing std MODTRAN atmos profiles
% emis_id2: PCRTM wavenumbers and stock emissivities for sensor id2

%% Output
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

%% Optimal estimation retrieval
converged = false;
if length(gammaFactor)< maxIter
    disp('Error: Not enough gamma factors for iteration!')
    ec = 1;
    return
end

m = length(y_obs);                                                        % length of measurement series
n = length(x0);                                                           % length of state series
K = zeros(m,n,maxIter);                                                   % list of Jacobians
x = zeros(n,maxIter+1);
y = zeros(m,maxIter+1);
A = zeros(n,n,maxIter);
d2 = zeros(maxIter,1);                                                    % convergence criteria
S_ep = zeros(m,m,maxIter);                                                % list of error covariance
S_aposteriori = zeros(n,n,maxIter);                                       % list of a posteriori covariance
dof = zeros(maxIter,1);                                                   % degree of freedom
chi2 = zeros(maxIter,2);                                                  % x-space; y-space; 
pvalue = zeros(maxIter,2);                                                % pvalue of chi-square test: x-space; y-space

%%
% deal with the first guess
x(:,1) = x0; 
y(:,1) = forwardPCRTM(pcrtmpath,ancpath,valid_rad,valid_emis,x0,lat,lon, ...
                      sfctype,pres,Ps,q,ozone,T,Ts,footprint_id, SRF, ...
                      MODTRAN_default_profile, emis_id2);

Sa_inv = inv(Sa);                                                         % inverted a priori covariance
% with no systematic bias b, 
S_epi = Sy;  
S_epi_inv = inv(S_epi);

for i = 1:maxIter
    xi = squeeze(x(:,i));
    Ki = JacobianPCRTM(pcrtmpath,ancpath,valid_rad,valid_emis,xi,lat,lon, ...
                       sfctype,pres,Ps,q,ozone,T,Ts,footprint_id, SRF, ...
                       MODTRAN_default_profile, emis_id2);
 
    % Using formulas from Turner and Lohnert 2014, doi:10.1175/JAMC-D-13-0126.1 
    B = gammaFactor(i)*Sa_inv + Ki'*S_epi_inv*Ki;                         % Equation (4)
    B_inv = inv(B);
    S_aposi = B_inv * (gammaFactor(i)^2*Sa_inv + Ki'*S_epi_inv*Ki) * B_inv;     % Equation (3)
    Gi = B_inv*Ki'*S_epi_inv;
    Ai = Gi*Ki;                                                           % Equation (5)
    dof(i) = trace(Ai);                                                   % Rodgers equation 2.80
    % Estimate next x and next y
    xii = xa + B_inv * Ki' * S_epi_inv * (y_obs-y(:,i)+Ki*(xi-xa));             % Equation (1)
    yii = forwardPCRTM(pcrtmpath,ancpath,valid_rad,valid_emis,xii,lat,lon, ...
                       sfctype,pres,Ps,q,ozone,T,Ts,footprint_id, SRF, ...
                       MODTRAN_default_profile, emis_id2);
    
    %%%%%%%%%%%%% Test convergence %%%%%%%%%%%%%%
    % cost function in x space Eq (5.28) in Rodgers 2000
       dx = x(:,i)-xii;       
       %chi2(i,1)= dx'* truncatedSVD(S_aposi) * dx;
       [chi2(i,1),pvalue(i,1)] = testChisquare(S_aposi,dx);
    % cost function in y space
       dy = y_obs-yii;
       % Eq 5.27 & 5.32 Rodgers 2000
       KSaKSep = Ki*Sa*Ki'+S_epi;
       S_deyd = S_epi*inv(KSaKSep)*S_epi;
       %chi2(i,2) = dy'* truncatedSVD(S_deyd) * dy;
       [chi2(i,2),pvalue(i,2)] = testChisquare(S_deyd,dy);
    % Test convergence   
    if n<= m   
       % in x-space
       % Convergence criterion Eq 5.29 Rodgers 2000
       d2i = chi2(i,1) ;
       d2i_limit = n./convergenceFactor;                                  
    else
        % in y-space
        % Convergence criterion Eq 5.33 Rodgers 2000
       dyy = y(:,i)-yii;
       d2i = dyy'*inv(S_deyd)*dyy;
       d2i_limit = m./convergenceFactor;
    end
    
    % Store the values
    x(:,i+1)=xii;
    y(:,i+1)=yii;
    K(:,:,i)=Ki;
    A(:,:,i)=Ai;
    d2(i) = d2i;
    S_ep(:,:,i)=S_epi;
    S_aposteriori(:,:,i)= S_aposi;
    
    if d2i<0
%        msg = 'WARNING: Negative convergence criterion';
%        fprintf('%s\n', msg);
        ec = 2;
        break
    end
    
    if i~=1
       if (abs(d2i)<=d2i_limit)&&(gammaFactor(i)==1)&&(d2i~=0)
           converged = true;
%           disp('Converged!')
           break
       elseif (i>2)&& (dof(i)==0)
           converged = false;
%           msg = 'WARNING: Zero degrees of freedom';
%           fprintf('%s\n', msg);
           ec = 3;
           break
%       else
%           disp('Not converged yet...')
       end
    end

end
%fprintf('Iterations performed: %d\n', i);

if converged
    iter = i;
    x_op = x(:,i);
    y_op = y(:,i);
    S_op = S_aposteriori(:,:,i);
    dgf = dof(i);
    K_op = squeeze(K(:,:,i));
else
    iter = i;
    x_op = x(:,i);
    y_op = y(:,i);
    S_op = S_aposteriori(:,:,i);
    dgf = dof(i);
    K_op = squeeze(K(:,:,i));
end

% calculate the nonlinearity parameter
% follow equations in Section 5.1 (Rodger 2000)
delta_x = sqrt(diag(S_op)) .* sign(y_op - y_obs);   
delta_y = y_op - y_obs - K_op*delta_x;
c2 = delta_y.* delta_y./ diag(Sy);
end
