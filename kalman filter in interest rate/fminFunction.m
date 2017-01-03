%======================================================================================
%***DESCRIPTION 
% This function calculated parameters in Kalman Filter to get prediction on interest rates.
% 
%***INPUTS
% data:  Price with different maturity dates
% R_est: Std of price with different maturity dates
% tau:   Maturity dates
%
%***OUTPUTS
% a(10*1):Parameters in Kalman Filter
%
%***CALL FUNCTION
% kalman_Log_Likelihood.m
%
%***REFERENCES 
% Forecasting, structural time series models and the Kalman filter ANDREW C. HARVEY, 2001 
% Linear Gaussian affine term structure models with unobservable factors:Calibration and yield forecasting, Date&Wang(2008)
% http://faculty.washington.edu/eeholmes/Files/Intro_to_kalman.pdf
%=========================================================================================

%%
%---------------------------------------------------------------------------------------
function [a] = fminFunction(data,R_est,tau)

y = data;

%Start with some reasonable initial parameter estimates
sigma_est_1_1 = R_est;
sigma_est_1_2 = R_est;
sigma_est_2_1 = R_est;
sigma_est_2_2 = R_est;

A_est_1 = 0.5;
A_est_2 = 0.5;

lambda_est_1 = 0.01;
lambda_est_2 = 0.01;

mu_est = 0;

% Set up a 10*1 matrix including all the parameters with initial values
startvals=[sigma_est_1_1;sigma_est_1_2;sigma_est_2_1;sigma_est_2_2;A_est_1;A_est_2;lambda_est_1;lambda_est_2;R_est;mu_est];

% "fminsearch" is a Nelder-Mead minimization matlab function
% Function options needed to set up in case that function would not stop 
options = optimset('fminsearch');
options = optimset(options,'Display','iter'); 

% Calculate the matrix "a" using "fminsearch"
% Call kalman_Log_Likelihood 
a = fminsearch('kalman_Log_Likelihood',startvals,options,y,[y(1),y(1)],[R_est,R_est],tau);
end
%---------------------------------------------------------------------------------------

