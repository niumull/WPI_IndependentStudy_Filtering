%%
% Description: 
% This file is to mimic the example of Auxiliary Particle Filter in 
%   A Tutorial on Particle Filters:
%       M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%            Section V-B.2
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
% w_k_j: Weights 
%
% Call functin
% resample.m
% f_x.m
%
% Notes
% The larger N we choose, the more explicit result we can get.
%
%%
function [x_k_j,w_k_j] = ASIR_PF(T,N,Q,R,z_k)
% %% 
% T=100;
% N= 100;
% Q= 1; 
% R=1; 
% x_k = zeros(1,T+1); % Initialized value
% Z_k = zeros(1,T+1);
% 
% for k = 2:T+1
%     x_k(k) = 0.5*x_k(k-1)+(25*x_k(k-1))/(1+x_k(k-1)^2)+8*cos(1.2*k)+sqrt(Q)*randn(1);
%     Z_k(k) = (x_k(k)^2)/20 + randn(1)*sqrt(R);
% end
% z_k = Z_k(2:T+1);

%%
% Initialize values
x_k_i = zeros(N,T+1); mu_k_i = zeros(N,T); mu_k = zeros(N,T); P_Z_mu_0 = zeros(N,T); 
P_Z_mu = zeros(N,T); P_Z = zeros(N,T); P_Z_x = zeros(N,T); w_k_i = ones(N,T+1)./N; 
x_k_j = zeros(N,T); w_k_j = zeros(N,T); H_k = zeros(N,T);

for k = 1:T
    
    %% Calculating the new step of x_k_i, s.t.x_k_j(:,k) = x_k_i(:,k+1)   
    % Sample mean
     mu_k_i(:,k) = f_x(x_k_i(:,k),k+1);
     rn = sqrt(Q)*randn(N,1);
     mu_k(:,k) = (mu_k_i(:,k).^2+rn.^2)/20;
    
    %% Calculate & normalize P(z_k|mu_k_i)
    var(:,k) = (mu_k_i(:,k).^4+6*(mu_k_i(:,k).^2)*Q+3*Q^2)/400-(mu_k(:,k).^2)+R;
    P_Z_mu_0(:,k) = exp((-(mu_k(:,k)-z_k(k)).^2)./(2*var(:,k))); % Calculate
    P_Z_mu(:,k) = P_Z_mu_0(:,k)./sum(P_Z_mu_0(:,k)); % Normalize
    
    %% Weights & Resample
    w_k_i(:,k+1) = w_k_i(:,k).*P_Z_mu(:,k);
    [~,~,i_j] = resample(x_k_i(:,k),w_k_i(:,k+1),N);
    
    for j = 1:N
        x_k_j(j,k) = f_x(x_k_i(i_j(j),k),k+1)+sqrt(Q)*randn(1); 
    end
    x_k_i(:,k+1) = x_k_j(:,k);
    
    %% Recalculating weights using SIR Filter
    H_k(:,k) =  (x_k_j(:,k).^2)/20;
    P_Z(:,k) = exp((-(H_k(:,k)-z_k(k)).^2)/(2*R));
    P_Z_x(:,k) = P_Z(:,k)./sum(P_Z(:,k));
    w_k_j(:,k) = P_Z_x(:,k);
end