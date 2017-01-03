function [X_new,F] = f_x(X_old,t)
%% 
% Description:
% f_x function is the model used in Unscented Filter
% 
% Input: value of predicted mean at time k-1, time step t
% 
% Output: value of predicted variance at time k
%
%% 
% Set value of parameters in function
R_s = 0.18; R_r = 0.15; M = 0.068; L_s = 0.0699; L_r = 0.0699; J = 0.0586;
T_l = 10; p = 1; h = 0.0001;
T_r = L_r/R_r;
sigma = 1 - M^2/(L_s*L_r);
K = M/(sigma*L_s*L_r);
gamma = R_s/(sigma*L_s) + (R_r*M^2)/(sigma*L_s*L_r^2);

% X_hat(:,1) = [200;200;50;50;350];
% X_hat_new = zeros(2*L+1,T); 
% Chi = zeros(5,2*L+1);
% interval = zeros(5,2*L+1);
% P_hat = 100^2.*eye(5);
u_1 = 350*cos(0.003*(t));
u_2 = 300*sin(0.003*(t));

% write function f_x
X_new(1) = X_old(1)+h*(-gamma*(X_old(1))+(K/T_r)*X_old(3)+K*p*X_old(5)*X_old(4)+(1/(sigma*L_s))*u_1);
X_new(2) = X_old(2)+h*(-gamma*(X_old(2))+(K/T_r)*X_old(4)-K*p*X_old(5)*X_old(3)+(1/(sigma*L_s))*u_2);
X_new(3) = X_old(3)+h*((M/T_r)*(X_old(2))-(1/T_r)*X_old(3)-p*X_old(5)*X_old(4));
X_new(4) = X_old(4)+h*((M/T_r)*(X_old(2))-(1/T_r)*X_old(4)+p*X_old(5)*X_old(3));
X_new(5) = X_old(5)+h*(((p*M)/(J*L_r))*(X_old(3)*X_old(2)-X_old(4)*X_old(1))-T_l/J);
F = [1-h*gamma,0,h*K/T_r,h*K*p*X_old(5),h*K*p*X_old(4);
     0,1-h*gamma,-h*K*p*X_old(5),h*K/T_r,-h*K*p*X_old(3);
     0,h*M/T_r,1-h/T_r,-h*p*X_old(5),-h*p*X_old(4);
     0,h*M/T_r,h*p*X_old(5),1-h/T_r,h*p*X_old(3);
     -h*p*M/(J*L_r)*X_old(4),h*p*M/J/L_r*X_old(3),h*p*M/J/L_r*X_old(2),-h*p*M/J/L_r*X_old(1),1];
end