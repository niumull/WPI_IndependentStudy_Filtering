%%
% Description: 
% This file is to calculate the RMSE of four kinds of Particle Filters to compare 
%   the prediction accuracy.
% 
% Reference:
%   A Tutorial on Particle Filters:
%       M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%
% Input:
% T: # of time steps
% N: # of iterations
% Q: Variance of estimator x
% R: Variance of observation Z
%
% Output:
% RMSE(Root Mean Square Error) of:
% SIR Particle filter
% Regularized Particle Filter
% Auxiliary Particle Filter
% Likelihood Particle Filter
%
% Call functin
% resample.m
% Regularized_PF.m
% SIR_PF.m
% ASIR_PF.m
% Likelihood_PF.m
%
% Notes
% The larger N we choose, the more explicit result we can get.
%
%%
function [Z_k_S,Z_k_R,Z_k_A,Z_k_L] = Compare(T,N,Q,R,z_k)

%% adjusted estimators 
[x_k_R,~] = Regularized_PF(T,N,Q,R,z_k);
[x_k_S,~] = SIR_PF(T,N,Q,R,z_k);
[x_k_A,~] = ASIR_PF(T,N,Q,R,z_k);
[x_k_L,~] = Likelihood_PF(T,N,Q,R,z_k);

% Observation at T+1 based on model we use
x_k_S_pre = mean(f_x(x_k_S(:,T),T+1)+sqrt(Q)*randn(N,1));
x_k_R_pre = mean(f_x(x_k_R(:,T+1),T+1)+sqrt(Q)*randn(N,1));
x_k_A_pre = mean(f_x(x_k_A(:,T),T+1)+sqrt(Q)*randn(N,1));
x_k_L_pre = mean(f_x(x_k_L(:,T),T+1)+sqrt(Q)*randn(N,1));
% x_k_S_pre = 0.5*x_k_S(:,T)+(25*x_k_S(:,T))./(1+x_k_S(:,T).^2)+8*cos(1.2*(T+1))+sqrt(Q)*randn(N,1);
% x_k_R_pre = 0.5*x_k_R(:,T)+(25*x_k_R(:,T))./(1+x_k_R(:,T).^2)+8*cos(1.2*(T+1))+sqrt(Q)*randn(N,1);
% x_k_A_pre = 0.5*x_k_A(:,T)+(25*x_k_A(:,T))./(1+x_k_A(:,T).^2)+8*cos(1.2*(T+1))+sqrt(Q)*randn(N,1);
% x_k_L_pre = 0.5*x_k_L(:,T)+(25*x_k_L(:,T))./(1+x_k_L(:,T).^2)+8*cos(1.2*(T+1))+sqrt(Q)*randn(N,1);

% Estimated Z_k at time T+1 
Z_k_S = (x_k_S_pre^2)/20 + randn(1)*sqrt(R);
Z_k_R = (x_k_R_pre^2)/20 + randn(1)*sqrt(R);
Z_k_A = (x_k_A_pre^2)/20 + randn(1)*sqrt(R);
Z_k_L = (x_k_L_pre^2)/20 + randn(1)*sqrt(R);