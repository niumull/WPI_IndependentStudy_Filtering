function [par] = fmin(Initials)
%% 
% Description:
% Estimate paremters in CIR and S&P model
%
% Input: initial values of 8 parameters
% Output: estimated values of 8 parameters
%
% Function: fminsearch

%% Set initial values
Q_sqrt = [Initials(1),Initials(2)];
k = Initials(3);
theta = Initials(4);
sigma = Initials(5);
mu = Initials(6);
a = Initials(7);
rho = Initials(8);

startvals = [Q_sqrt(1),Q_sqrt(2),k,theta,sigma,mu,a,rho];

% Do fminsearch
options = optimset('fminsearch');
options = optimset(options,'Display','iter','MaxIter',1000); 
par = fminsearch('MLL_fmin',startvals,options);
par(1)=par(1)^2;
par(2)=par(2)^2;
end
