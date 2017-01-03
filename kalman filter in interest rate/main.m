%======================================================================================
%***DESCRIPTION
% This file is the main file to predict interest rate. 
%  (1) Collect data of 10-Years, 8-Years, 4-Years, 2-Years, 1-Years UK bond yields 
%  (2) Use fminFunction.m file to calculate parameters 
%  (3) Calculate weight of price with these 5 different maturity dates.
%  (4) Use estimation.m file to calculate estimated interest rate with different maturity dates seperately
%  (5) Add weight in the estimation
% The result is represented by y_hat
% 
%***CALL FUNCTION
% fminFunction.m
% estimation.m
%
%***INPUTS
% data_@ represents the bond yield with maturity date @ years
% R_est_@ represents the standard deviation of a bond yield with maturity date @ years
%=========================================================================================

%---------------------------------------------------------------------------------------
% Input price of different maturity dates in data_10,data_8,data_4,data_2,data_1                     
% Input std of price with different maturity dates in R_est_10,R_est_8,R_est_4,R_est_2,R_est_1
data_10 = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','B3:B235');
R_est_10 = sqrt(xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','B237'));
data_8 = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','C3:C235');
R_est_8 = sqrt(xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','C237'));
data_4 = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','D3:D235');
R_est_4 = sqrt(xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','D237'));
data_2 = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','E3:E235');
R_est_2 = sqrt(xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','E237'));
data_1 = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','F3:F235');
R_est_1 = sqrt(xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Price','F237'));
%---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
% Calculated parameters with different maturity dates
% Call fminFunction
[v_10] = fminFunction(data_10,R_est_10,10);
[v_8]  = fminFunction(data_8,R_est_8,8);
[v_4]  = fminFunction(data_4,R_est_4,4);
[v_2]  = fminFunction(data_2,R_est_2,2);
[v_1]  = fminFunction(data_1,R_est_1,1);
%---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
% Calculated weights based on contribution of variance
% Using the first five columns
UKChange = xlsread('E:\16 Summer Research\Linear Gaussian.xlsx','Change','B3:F235');
% Covariance matrix
UK_Cov = cov(UKChange);
UK_Corr = corrcov(UK_Cov);
% Eigen value
UK_Eig = eig(UK_Corr);
% Weights
UK_weights = UK_Eig/sum(UK_Eig);
%---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
% Calculated estimated bond yields with different maturity dates
% Call estimation
y_hat_10 = estimation(v_10,data_10,10);
y_hat_8 = estimation(v_8,data_8,8);
y_hat_4 = estimation(v_4,data_4,4);
y_hat_2 = estimation(v_2,data_2,2);
y_hat_1 = estimation(v_1,data_1,1);
%---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
% Weighted bond yield y_hat which represented estimated interest rate
y_hat = y_hat_10*UK_weights(1)+y_hat_8*UK_weights(2)+y_hat_4*UK_weights(3)+y_hat_2*UK_weights(4)+y_hat_1*UK_weights(5);
%---------------------------------------------------------------------------------------
