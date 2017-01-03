clear all
clc
%% Description:
% Calculate estimated parameters in CIR and S&P model
% presented in par (1*8)
% do test work using two methods

%% Initials
rng(100,'twister')
% Input initial values
% Initials = [Q1,Q2,k,theta,sigma,mu,a,rho]
Initials = zeros(1,8);
Initials(1) = 0.5^(0.25); %Q1=sqrt(0.5), to keep Q1 positive
Initials(2) = 0.5^(0.25);
Initials(3) = sqrt(0.088);%0.088
Initials(4) = sqrt(0.035);%0.035
Initials(5) = 1;
Initials(6) = 1.66;
Initials(7) = 0.031;
Initials(8) = -0.84;

%% observation data
Y_obs(:,1) = xlsread('Shiller data.xlsx','Sheet3','F3:F66');
Y_obs(:,2) = xlsread('Shiller data.xlsx','Sheet3','E3:E66');

%% Put initial values into ParaEstimate.m function, estimate parameters
%  par_fmin = fmin(Initials);

%% Put initial values into newton.m function, estimate parameters
tol = 0.001; nmax = 100;h = 2;

par_newton = newton(Initials,tol,nmax,h,Y_obs); %new parameters show Q1 & Q2 instead of sqrt(Q)

%% Use initial values to estimate 
% par_initial = Initials;

%% Test results
Y_real(:,1) = xlsread('Shiller data.xlsx','Sheet3','F67:F72');
Y_real(:,2) = xlsread('Shiller data.xlsx','Sheet3','E67:E72');

% [Y_estimate_09,X_estimate_09,rmse_09] = test_09(par_fmin,Y_real);
% [Y_estimate_rolling,X_estimate_rolling,rmse_rolling] = test_rolling(par_fmin,Y_real);
% 
[Y_estimate_09,X_estimate_09,rmse_09] = test_09(par_newton,Y_real);
[Y_estimate_rolling,X_estimate_rolling,rmse_rolling] = test_rolling(par_newton,Y_real);

%  [Y_estimate_09,X_estimate_09,rmse_09] = test_09(par_initial,Y_real);
%  [Y_estimate_rolling,X_estimate_rolling,rmse_rolling] = test_rolling(par_initial,Y_real);

figure(1)
hold on
 plot([2010:1:2015],exp(Y_estimate_rolling(:,1)));
 plot([2010:1:2015],exp(Y_real(:,1))); 
%
%ylim([-1 1])
legend('Estimated Value','Observations');
title('Comparison Between Estimated Dividend Yield and Observations')
xlabel('Year')
ylabel('Dividend Yield')
grid ;
hold off

% 
figure(2)
hold on
plot([2010:1:2015],Y_estimate_rolling(:,2));
plot([2010:1:2015],Y_real(:,2));
%ylim([-1 1])
legend('Estimated Value','Observations');
title('Comparison Between Estimated S&P Returns and Observations')
xlabel('Year')
ylabel('S&P Returns')
grid ;
hold off

