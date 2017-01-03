function result = runParticleSimulationSmoother(mStates, mWeights, vY)

mLogVols1 = mStates(2,:,:);
mLogVols2 = mStates(3,:,:);
mAlpha = mStates(1,:,:);
iNumObs = size(vY,2);
	
vWeights = mWeights(:,iNumObs)';
index = randm(vWeights,1);
dLogVols1 = mLogVols1(1,index,iNumObs);
dLogVols2 = mLogVols2(1,index,iNumObs);
dAlpha = mAlpha(1,index,iNumObs);

% ===================================================== //
% store the draw
% ===================================================== //
dSignalToNoise = exp(dLogVols2)./exp(dLogVols1);
dTheta = (sqrt(dSignalToNoise.^2 + 4*dSignalToNoise) - 2 - dSignalToNoise)/2;
mSmoothed = zeros(4,iNumObs);
mSmoothed(1,iNumObs) = dAlpha;
mSmoothed(2,iNumObs) = exp(dLogVols1/2);
mSmoothed(3,iNumObs) = exp(dLogVols2/2);
mSmoothed(4,iNumObs) = dTheta;

for t=(iNumObs-1):-1:1
	% ===================================================== %
	% compute backwards weights
	% ===================================================== %
	vLogWeights = logSmoothedWeights(dLogVols1, dLogVols2, dAlpha, mLogVols1(1,:,t), mLogVols2(1,:,t), mAlpha(1,:,t));

	% ============================================== %
    % Normalize the importance weights
   	% ============================================== %
	vLogWeights = vLogWeights + log(mWeights(:,t)');
	vWeights = exp(vLogWeights - max(vLogWeights));
    vWeights = vWeights/sum(vWeights);
		
	% ===================================================== %
	% draw a candidate backwards
	% ===================================================== %
	index = randm(vWeights,1);
	dLogVols1 = mLogVols1(1,index,t);
	dLogVols2 = mLogVols2(1,index,t);
	dAlpha = mAlpha(1,index,t);

	% ===================================================== //
	% store the draw
	% ===================================================== //
	dSignalToNoise = exp(dLogVols2)./exp(dLogVols1);
	dTheta = (sqrt(dSignalToNoise.^2 + 4*dSignalToNoise) - 2 - dSignalToNoise)/2;
    mSmoothed(1,t) = dAlpha;
	mSmoothed(2,t) = exp(dLogVols1/2);
	mSmoothed(3,t) = exp(dLogVols2/2);
	mSmoothed(4,t) = dTheta; 
end
result = mSmoothed;
return

function vLogWeights = logSmoothedWeights(dLogVols1,  dLogVols2, dAlpha, mLogVols1, mLogVols2, mAlpha)
	
vLogWeights = 0;
% first distribution
vLogWeights = vLogWeights  - 0.5*log(2*pi) - 0.5*log(0.2) - 0.5*(1/0.2)*(dLogVols1 - mLogVols1).^2;
vLogWeights = vLogWeights  - 0.5*log(2*pi) - 0.5*log(0.2) - 0.5*(1/0.2)*(dLogVols2 - mLogVols2).^2;
vLogWeights = vLogWeights - 0.5*log(2*pi) - 0.5*log(exp(mLogVols2)) -0.5*(1 ./exp(mLogVols2)).*(dAlpha - mAlpha).^2;
return

function result = randm(probs,m)
% ======================================== %
% multinomial random number generator
% ======================================== %
result = zeros(m,1);
for j=1:m
    result(j,1) = length(find(cumsum(probs')<rand))+1; 
end
return