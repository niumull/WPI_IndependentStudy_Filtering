
clear all

MC = 300;
T = 100;
N = 50;
Q = 10; 
R = 1; 

[RMSE_SIR_MC,RMSE_RPF_MC,RMSE_ASIR_MC,RMSE_LPF_MC] = Main(MC,T,N,Q,R);
% Calculating RMSE for MC times