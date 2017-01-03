function [Y_hat_pre] = ukf_modified(H_k,Q,R,T,L,W0,I,Y)
%%
% Description: 
% To accomplish the prediction of non-linear system:
% X_k = f_x(X_(k-1))+w_k
% Y_k = H_k*X_k+v_k
% where w_k~N(0,Q),v_k~N(0,R)

% Set initial value of X, mu and sigma to produce sigma points
d1 = length(H_k); 
d2 = 2;
X_hat(:,1) = [200;200;50;50;350];% mu
P_hat = 100^2.*eye(d1);% sigma
a = sqrt(1/(1-W0));

% Calculate w by using a
for i = 1:2*L+1;
    if i == L+1
        w(i) = 1-1/a^2;
    else
        w(i) = 1/(2*L*a^2);
    end
end

Q_k = Q*eye(d1); 
R_k = R*eye(d2);
delta_Q_k = I*eye(d1);% input delta_Q_k to improve accuracy
for t = 1:T+1
    [X_hat(:,t+1),P_hat] = update(X_hat(:,t),P_hat,Y(:,t),L,a,w,Q_k,delta_Q_k,H_k,R_k,t);
end

Y_hat = X_hat(1:2,2:end)+sqrt(R)*randn(d2,T+1);
Y_hat_pre = Y_hat(:,end);
end





