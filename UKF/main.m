clear all;
clc
%% 
% Description:
% To accomplish the prediction of non-linear system:
% X_k = f_x(X_(k-1))+w_k
% Y_k = H_k*X_k+v_k
% where w_k~N(0,Q),v_k~N(0,R)

%% 
% Set initial values in function
H_k = [1,0,0,0,0;0,1,0,0,0];
d1 = length(H_k); 
d2 = 2;% length of Y
Q = 0.01^2; % X
R = 0.1^2; % Y
% Set T as time steps, set L such that sigma points' # is (2*L+1)  
T = 1800;
L = 50;
% Set W0 to produce sigma points
W0 = 1/3; 
I = 0.05^2;
MC = 2;
X(:,1) = zeros(d1,1);

% Set observations Y(d2*T+1), unsing function f_x.m
for t = 1:T+1
    X(:,t+1) = f_x(X(:,t),t)'+ sqrt(Q)*randn(d1,1);
end
X = X(:,2:end);% d1*T+1
Y = H_k*X+sqrt(R)*randn(d2,T+1);% d1*T+1
Y_obs = Y(:,end);%take it to calcuate rmse

% Calculate predicted Y(Y_hat_pre), using function ukf_modified.m
for i = 1:MC
     Y_hat_pre(:,i) = ukf_modified(H_k,Q,R,T,L,W0,I,Y);
end

% Calculate rmse(d2*1)
for k = 1:d2
rmse(k,:) = sqrt(sum((Y_hat_pre(k,:)-Y_obs(k)).^2./MC));
end

