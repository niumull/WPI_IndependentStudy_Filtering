%%
% Description: 
% This file is to replicate the calculation of the RMSE of four kinds of Particle 
%   Filters using Monte Carlo method to compare the prediction accuracy.
% 
% Reference:
% A Tutorial on Particle Filters:
%   M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%
% Input:
% MC: # of Monte Carlo runs
% T: # of time steps
% N: # of iterations
% Q: Variance of estimator x
% R: Variance of observation Z
%
% Output:
% Monte Carlo mean of RMSE(Root Mean Square Error) of:
% SIR Particle filter
% Regularized Particle Filter
% Auxiliary Particle Filter
% Likelihood Particle Filter
%
% Call functin
% resample.m
% f_x.m
% Regularized_PF.m
% SIR_PF.m
% ASIR_PF.m
% Likelihood_PF.m
% Compare.m
%
%%
function[RMSE_SIR_MC,RMSE_RPF_MC,RMSE_ASIR_MC,RMSE_LPF_MC] = Main(MC,T,N,Q,R)
% Calculating RMSE for MC times

rng(0,'twister'); 

%randn(100,1);
%% To ensure the accuracy, we use the same set of observations for each
% Initialized value
x_k = zeros(1,T+1); Z_k = zeros(1,T+1); Z_k_S = zeros(1,MC);
Z_k_R = zeros(1,MC); Z_k_A = zeros(1,MC); Z_k_L = zeros(1,MC);

for k = 2:T+2
    x_k(k) = f_x(x_k(k-1),k)+sqrt(Q)*randn(1);
    Z_k(k) = (x_k(k)^2)/20 + randn(1)*sqrt(R);
end
z_k = Z_k(2:T+1);
z_k_pre = Z_k(T+2);

for i = 1:MC
    [Z_k_S(i),~,~,~] = Compare(T,N,Q,R,z_k);    
end
RMSE_SIR_MC = sqrt(sum((Z_k_S-z_k_pre).^2)/MC);
disp('SIR Done.')

for i=1:MC
    [~,Z_k_R(i),~,~] = Compare(T,N,Q,R,z_k);
end
RMSE_RPF_MC = sqrt(sum((Z_k_R-z_k_pre).^2)/MC);
disp('RPF Done.')

for i=1:MC
    [~,~,Z_k_A(i),~] = Compare(T,N,Q,R,z_k);
end
RMSE_ASIR_MC = sqrt(sum((Z_k_A-z_k_pre).^2)/MC);
disp('ASIR Done.')

for i=1:MC
    [~,~,~,Z_k_L(i)] = Compare(T,N,Q,R,z_k);
end
RMSE_LPF_MC = sqrt(sum((Z_k_L-z_k_pre).^2)/MC);
disp('LPF Done.')

   