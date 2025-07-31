function [chi2,pvalue] = testChisquare(S,z)
% run the chi sqaure test following Rodgers 2000
% Chapter 12: Testing and Validating
%%% input %%%
% S: a real symmetric matrix
% z: a random vector which is assumed to have zero mean
%%% output %%%
% chi2: observed chi-square value
% pvalue: correponding p-value given observation and dof
% S = rand(7,7);         % for test
% z = rand(7,1);
%% 
[~,D,W] = eig(S);
lam = diag(D);
z_prime = W' * z;

% Deal with possible singular matrices
% Rodgers Chapter 12 Section 2
idx_notnull  = ( abs(lam) > 1e-4*max(abs(lam)) ) ;
% idx_notnull = (abs(lam) > 1e-5 );
dof = sum(idx_notnull);                                % degree of freedom: number of independent pieces of information
chi2_sequ = z_prime(idx_notnull).^2 ./ lam(idx_notnull); 
chi2 = sum(chi2_sequ);                                 % Eq.(12.1) in Rodgers 2000

% given observed chi2 and dof, determine p value
% p-value measures the probability that an observed differene could have
% occurred just by random chance
pvalue = 1 - chi2cdf(chi2, dof);
