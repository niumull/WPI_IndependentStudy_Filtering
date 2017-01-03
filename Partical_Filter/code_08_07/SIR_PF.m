 function [x_k_j,w_k_i] = SIR_PF(T,N,Q,R,z_k)

 % function [x_k_j,w_k_i] = SIR_PF(T,N,Q,R,z_k)
%%
% Description: 
% This file is to mimic the example of SIR Particle Filter in 
%   A Tutorial on Particle Filters:
%       M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%            Section V-B.1
%
% Input:
% T: # of time steps
% N: # of iterations
% Q: Variance of estimator x
% R: Variance of observation Z
% z_k: Observations
%
% Output:
% x_k_j: estimated samples
% w_k_i: Weights 
%
% Call functin
% resample.m
% f_x.m
%
% Notes
% The larger N we choose, the more explicit result we can get.
%


% Initialize values
x_k_i = zeros(N,T+1); H_k = zeros(N,T); P_Z = zeros(N,T); w_k_i = ones(N,T)./N;
x_k_j = zeros(N,T); w_k_j = zeros(N,T);


for k = 1:T
    
    % Calculating the new step of x_k_i, s.t.x_k_j(:,k) = x_k_i(:,k+1)
    x_k_j(:,k) = f_x(x_k_i(:,k),k+1)+sqrt(Q)*randn(N,1);
    % mean of x_k_j
    H_k(:,k) =  (x_k_j(:,k).^2)/20;
    
    % Normally distributed density
    % Canceled "(1/sqrt(2*pi*R))" as a constant
    P_Z(:,k) = exp((-(H_k(:,k)-z_k(k)).^2)/(2*R));
    % Normalize weights
    w_k_i(:,k) = P_Z(:,k)./sum(P_Z(:,k));
    
    % Resample estimators and weights to avoid degeneracy
    [x_k_j(:,k),w_k_j(:,k),~] = resample(x_k_j(:,k),w_k_i(:,k),N);
    % Overwright the new x_k_i & w_k_i
    x_k_i(:,k+1) = x_k_j(:,k);
    w_k_i(:,k) = w_k_j(:,k);
end
end



