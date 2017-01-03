%%
% Description: 
% This file is to mimic the example of 'Likelihood' Particle Filter in 
%   A Tutorial on Particle Filters:
%       M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%            Section A
%
% Input:
% T: # of time steps
% N: # of iterations
% Q: Variance of estimator x
% R: Variance of observation Z
% z_k: Observations
%
% Output:
% x_k_i: estimated samples
% w_k_i: Weights 
%
% Call functin
% resample.m
% f_x.m
%
% Notes
% The larger N we choose, the more explicit result we can get.
%
%%
function [x_k_i,w_k_i] = Likelihood_PF(T,N,Q,R,z_k)
% Initialize values
x_k_i = zeros(N,T); w_k_i = ones(N,T)./N; s_k_i = zeros(N,T); w_k = zeros(N,T);
mu_k_i = zeros(N,T); P_x = zeros(N,T); x_k_j = zeros(N,T); w_k_j = zeros(N,T); 
N_eff = zeros(1,T);

% Calculated x_k in time 1 
x_k_i(:,1) = 8*cos(1.2*1)+sqrt(Q)*randn(N,1);
for k = 2:T
    % s_k_i is (x_k_i)^2 in reference
    s_k_i(:,k) = (z_k(k)*ones(N,1)-randn(N,1)*sqrt(R))*20;
    % If received a nagative s_k, repeat the process 
    for i = 1:N
        while s_k_i(i,k) < 0
            s_k_i(i,k) = (z_k(k)-randn(1)*sqrt(R))*20;
        end
        % randomly set x_k to be half positive & half negative
        U = rand(N,1);
        if U(i) > 0.5
            x_k_i(i,k) = sqrt(s_k_i(i,k));
        else
            x_k_i(i,k) = -sqrt(s_k_i(i,k));
        end
    end
    % Calculate P(x_k_i|x_k-1_i) using x_k_i above and x_k-1_i from last time step 
    mu_k_i(:,k) = f_x(x_k_i(:,k-1),k);
    P_x(:,k) = exp((-(mu_k_i(:,k)-x_k_i(:,k)).^2)/(2*Q));
    % Weights
    w_k(:,k) = w_k_i(:,k-1).*P_x(:,k).*x_k_i(:,k);
    w_k_i(:,k) = w_k(:,k)./sum(w_k(:,k));
    % Checking degeneracy problem, if exist, resample
    N_eff(k) = 1/(sum(w_k_i(:,k).^2));
    if N_eff < N
        [x_k_j(:,k),w_k_j(:,k),~] = resample(x_k_i(:,k),w_k_i(:,k),N);
    end
    % Overwrite estimators and weights
    x_k_i(:,k) = x_k_j(:,k);
    w_k_i(:,k) = w_k_j(:,k);
end

