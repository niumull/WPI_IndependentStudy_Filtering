
clc;
clear all;
close all;

% Fix the see on the random number generator
randn('state',123);
rand('state',123);

iN = 10000; % Number of particles
iT = 1000; % Number of observations 

mu = 0.5;
phi = 0.975;
sigma2 = 0.04;
	
vParms = [phi ; sigma2 ; mu ];

vEta = randn(1,iT+100);
vEps = randn(1,iT);
   
vLogVolsTrue = zeros(1,iT+100);

for t=1:(iT+100)
    if(t==1)
        vLogVolsTrue(t) = randn(1,1);
    else
		vLogVolsTrue(t) = mu + phi*(vLogVolsTrue(t-1) - mu) + sqrt(sigma2)*vEta(t);
    end
end

% NOTE: This code calls a MEX file for the resampling step.
% The MEX file has been compiled for 64 bit machines
% If you get an error, it is because the MEX file will not run 
% on your computer.
% To run the file, change the Matlab code on LINE 61 of
% raoBlackwellizedAPF.m to run the resampling algorithm written in Matlab.

vLogVolsTrue = vLogVolsTrue(:,101:end);
   
vY = exp(vLogVolsTrue/2).*vEps;

result = sisrFilterStochasticVolatility(vY, vParms, iN);

vLogVolsSISR = result.logFilter;
dLoglikeSISR = result.loglike;

result = bootstrapFilterStochasticVolatility(vY, vParms, iN);

vLogVolsBF = result.logFilter;
dLoglikeBF = result.loglike;

result = auxiliaryParticleFilterStochasticVolatility(vY, vParms, iN);

vLogVolsAPF = result.logFilter;
dLoglikeAPF = result.loglike;

figure(1);
plot(1:iT,vLogVolsTrue,'r',1:iT,vLogVolsBF,'b',1:iT,vLogVolsSISR,'g',1:iT,vLogVolsAPF,'y');
legend('True state variable','bootstrap filter','SISR filter','Auxiliary particle filter');
%plot(1:iT,vLogVolsTrue,'r',1:iT,vLogVolsBF,'b',1:iT,vLogVolsSISR,'g');