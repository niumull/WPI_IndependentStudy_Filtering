clc;
clear all;
close all;

% Fix the see on the random number generator
randn('state',123);
rand('state',123);

vInflation = load('inflation.mat');

vY = vInflation.inflation';
iN = 10000; % Number of particles
[iNumDim, iT] = size(vY); % Number of observations 

	
sigma2_trans = 0.2; % Values chosen by Stock and Watson (2007)
sigma2_meas = 0.2;

vParms = [sigma2_meas ; sigma2_trans];

% NOTE: This code calls a MEX file for the resampling step.
% The MEX file has been compiled for 64 bit machines
% If you get an error, it is because the MEX file will not run 
% on your computer.
% To run the file, change the Matlab code on LINE 61 of
% raoBlackwellizedAPF.m to run the resampling algorithm written in Matlab.
        
result = raoBlackwellizedAPF(vY, vParms, iN);

mFilter = result.filter;
mPredictor = result.predictor;
dLoglike = result.loglike;
mStates = result.states;
mWeights = result.weights;

figure(1)
plot(1:iT,vY,'r',1:iT,mFilter(1,:),'b');
legend('Observed inflation','Filtered trend');

figure(2)
plot(1:iT,vY,'r',1:iT,mPredictor(1,:),'b');
legend('Observed inflation','One-step ahead predicted trend');

figure(3)
plot(1:iT,mFilter(2,:),'b');
legend('Filtered volatility of the irregular component');

figure(4)
plot(1:iT,mFilter(3,:),'b');
legend('Filtered volatility of the trend');

figure(5)
plot(1:iT,mFilter(4,:),'b');
legend('Filtered value of the MA(1) parameter');

figure(6)
plot(1:iT,mFilter(5,:),'b');
legend('Filtered value of the signal to noise ratio');

return;

iNumDraws = 1000;
mSmoothed = zeros(4,iT);

% Smoothing which is slower to run
for j=1:iNumDraws
    result = runParticleSimulationSmoother(mStates, mWeights, vY);
    mSmoothed = mSmoothed + result/iNumDraws;
end

figure(7)
plot(1:iT,mPredictor(2,:),'g',1:iT,mFilter(2,:),'r',1:iT,mSmoothed(2,:),'b');

figure(8)
plot(1:iT,mPredictor(3,:),'g',1:iT,mFilter(3,:),'r',1:iT,mSmoothed(3,:),'b');

figure(9)
plot(1:iT,mPredictor(4,:),'g',1:iT,mFilter(4,:),'r',1:iT,mSmoothed(4,:),'b');



