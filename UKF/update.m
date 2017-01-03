function [X_new,P_new] = update(X_old,P_old,Y,L,a,w,Q_k,delta_Q_k,H_k,R_k,t)
%%
% Description:
% using ukf to update predicted x
% Input:
% X_old,P_old: value of predicted mean & variance at time t-1 
% Y:observations 
% L:time step
% a:tuning paremeter to spread sigma points, w is a function corresponding to a
% Q_k,R_k, H_k:parameters in the nonlinear system
% delta_Q_k: improve accuracy
% t: time step
% Output:
% X_new, P_new: updated predicted mean and variance
%% 
% Calculate predicted mean(X_hat) and variance by using sigma points
d1 = length(H_k); 
d2 = 2;
Chi = zeros(d1,2*L+1);
interval = zeros(d1,2*L+1);
X_hat = zeros(d1,1); 
for j = 1:d1
    interval(j,:) = linspace((X_old(j)-a*sqrt(L*P_old(j))),(X_old(j)+a*sqrt(L*P_old(j))),2*L+1);
    for i = 1:2*L+1;
        [Chi(:,i),~] = f_x(interval(:,i),i);
    end
    X_hat(j) = sum(w.*Chi(j,:));
end

[~,F_k] = f_x(X_old,t);

% Update by using modified unscented filter
telta_x_k = X_old-X_hat;
DELTA_P = F_k*(telta_x_k*telta_x_k'+Q_k)*F_k'-F_k*P_old*F_k';
P_K_k = F_k*P_old*F_k'+DELTA_P+Q_k;

sample = 0;
for i = 1:2*L+1;
    sample = sample+w(i)*(Chi(:,i)-X_hat)*(Chi(:,i)-X_hat)';
end

delta_P = P_K_k-(sample+Q_k);
Q_hat_k = DELTA_P+Q_k+delta_P+delta_Q_k;
P_hat_K_k = F_k*P_old*F_k'+Q_hat_k;
 
y_hat = [X_hat(1);X_hat(2)];
P_y_y = H_k*P_hat_K_k*H_k'+R_k;
P_x_y = P_hat_K_k*H_k';
X_new = X_hat+P_x_y*(P_y_y^(-1))*(Y-y_hat);
P_new = abs(P_hat_K_k-P_x_y*(P_y_y^(-1))*P_x_y');
end





