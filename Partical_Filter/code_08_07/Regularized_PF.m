%%
% Description: 
% This file is to mimic the example of Regularised Particle Filter in 
%   A Tutorial on Particle Filters:
%       M. Sanjeev Arulampalam, Simon Maskell, Neil Gordon, and Tim Clapp
%            Section V-B.3
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
%%
function [x_k_i,w_k_i] = Regularized_PF(T,N,Q,R,z_k)
% Initialize values
x_k_i = zeros(N,T+1); H_k = zeros(N,T); P_Z = zeros(N,T); w_k_i = ones(N,T)./N; 
x_k_j = zeros(N,T); w_k_j = zeros(N,T); N_eff = zeros(1,T); x_Reg = ones(N,T); 
xregM = zeros(1,T); xregm = zeros(1,T); x_norm = zeros(N,T); epsilon = zeros(1,T); 

% Constants in equation (77) (78)
% C_n_x = (pi^(N/2))/factorial(N/2);
% A = ((8/C_n_x)*(N+4)*((2*sqrt(pi))^N))^(1/(N+4));
A = 4/(N+2)^(1/(N+4));
h_opt = A*N^(-1/(N+4));

% Regularized Partical Filter
for k = 1:T
    % Calculating the new step of x_k_i, s.t.x_k_j(:,k) = x_k_i(:,k+1)
    x_k_j(:,k) = f_x(x_k_i(:,k),k+1)+sqrt(Q)*randn(N,1);
    H_k(:,k) =  (x_k_j(:,k).^2)/20;
    P_Z(:,k) = exp((-(H_k(:,k)-z_k(k)).^2)/(2*R));
    w_k_i(:,k) = P_Z(:,k)./sum(P_Z(:,k));
    % Empirical covariance matrix of x_k_i
    S_K = cov(x_k_j(:,k));
    if S_K > 0
        D_k = chol(S_K); 
    else
        D_k = 0;
    end
    N_eff(k) = 1/(sum(w_k_i(:,k).^2));
    if N_eff(k) < N
        [x_k_j(:,k),w_k_j(:,k),~] = resample(x_k_j(:,k),w_k_i(:,k),N);
        % Draw epsilon from the Epanechnikov Kernel (76)
        xregM(k) = max(x_k_j(:,k))-std(x_k_j(:,k));
        xregm(k) = min(x_k_j(:,k))+std(x_k_j(:,k)); 
        x_Reg(:,k) = linspace(xregM(k),xregm(k),N);
        for i = 1:N
            x_norm(i,k) = norm(1/(h_opt*det(D_k))*(x_Reg(i,k)-x_k_j(i,k)));  
            if x_norm(i,k) < 1  
%                 epsilon(i,k) = ((N+2)/(2*C_n_x))*(1-x_norm(i,k)^2);%0<epsilon<1?   
                  epsilon(i,k) =  exp(-(x_Reg(i,k)-x_k_j(i,k))'*(x_Reg(i,k)-x_k_j(i,k)))/(2*D_k^2);
            else
                epsilon(i,k) = 0;
            end
        end   
        % Get new x_k_i
        x_k_i(:,k+1) = x_k_j(:,k)+h_opt*D_k*epsilon(:,k);   
        w_k_i(:,k) = w_k_j(:,k);
    end
end
end



