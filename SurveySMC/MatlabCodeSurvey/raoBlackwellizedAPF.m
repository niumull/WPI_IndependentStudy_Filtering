function result = raoBlackwellizedAPF(vY, vParms, iN)

[iDim,iNumObs] = size(vY);

% ==============================
% Storage Matrices
% ==============================
mPredictedEstimates = zeros(5,iNumObs);
mFilteredEstimates = zeros(5,iNumObs);
mWeights = zeros(iN,iNumObs);
mStates = zeros(3,iN,iNumObs);

vESS = zeros(1,iNumObs);

% ============================================== //
% Initialize the particles
% ============================================== //
[mLogVols, mFilter] = initializeParticles(iN);
dLoglike = 0;
vWeights = ones(1,iN)*(1/iN);

for t=1:iNumObs
		% ============================================== //
    	% Run the prediction step of the Kalman Filter
    	% ============================================== //
		vSignalToNoise = exp(mLogVols(2,:))./exp(mLogVols(1,:));
		vTheta = (sqrt(vSignalToNoise.^2 + 4*vSignalToNoise) - 2 - vSignalToNoise)/2;
		mPredictions = KalmanPredictionStep(mLogVols, mFilter, iN);
			  
		mPredictedEstimates(1,t) = sum(vWeights.*mPredictions(1,:));
		mPredictedEstimates(2,t) = sum(vWeights.*exp(mLogVols(1,:)/2));
	    mPredictedEstimates(3,t) = sum(vWeights.*exp(mLogVols(2,:)/2));
		mPredictedEstimates(4,t) = sum(vWeights.*vTheta);
        mPredictedEstimates(5,t) = sum(vWeights.*vSignalToNoise);

		% ============================================== //
    	% Compute the importance weights
    	% ============================================== //
		vLogWeights = logIncrementalImportanceWeights(mLogVols, mPredictions, vY(t));

		% ============================================== %
    	% Normalize the importance weights
   		% ============================================== %
		dMaxWeight = max(vLogWeights);
		vWeights = exp(vLogWeights - dMaxWeight);
		dLoglike = dLoglike + log(mean(vWeights)) + dMaxWeight;
    	vWeights = vWeights/sum(vWeights);
	   
		% ============================================== % 
		% Compute the ESS
		% ============================================== %
		vESS(1,t) = 1/sum(vWeights.^2);
		
		% ============================================== %
    	% Resample
   		% ============================================== %
        % NOTE: If you get an error here it is likely that the MEX file which contains the resampling step is not working
        % on your computer. You can always use the first resampling algorithm
        % which is coded in Matlab
        
		%vIndex = systematicR_matlab(1:iN,vWeights');
        
        %vIndex = systematicResampling(vWeights);
        %vIndex = multinomialResampling(vWeights); 
        %vIndex = stratifiedResampling(vWeights);
        vIndex = residualResampling(vWeights);
        
		mLogVols = mLogVols(:,vIndex);
		mFilter = mFilter(:,vIndex);
		vWeights = ones(1,iN)*(1/iN);

		% ============================================== %
        % Draw new log volatilities
    	% ============================================== %
		mLogVols = drawVolatilities(vParms, mLogVols, iN);

		% ============================================== %
    	% Run the Kalman prediction step again
   		% ============================================== %
		mPredictions = KalmanPredictionStep(mLogVols, mFilter, iN);
		
		% ============================================== %
    	% Run the update step of the Kalman filter
   		% ============================================== %
		mFilter = KalmanUpdateStep(mLogVols, mPredictions, iN, vY(t));

		% ============================================== %
    	% Estimate the state variables
   		% ============================================== %
		vSignalToNoise = exp(mLogVols(2,:))./exp(mLogVols(1,:));
		vTheta = (sqrt(vSignalToNoise.^2 + 4*vSignalToNoise) - 2 - vSignalToNoise)/2;
		mFilteredEstimates(1,t) = sum(vWeights.*mFilter(1,:));
		mFilteredEstimates(2,t) = sum(vWeights.*exp(mLogVols(1,:)/2));
		mFilteredEstimates(3,t) = sum(vWeights.*exp(mLogVols(2,:)/2));
        mFilteredEstimates(4,t) = sum(vWeights.*vTheta);
        mFilteredEstimates(5,t) = sum(vWeights.*vSignalToNoise);
		
        % ============================================== %
    	% Store the particles for smoothing later
   		% ============================================== %
		mStates(1,:,t) = mFilter(1,:)';
		mStates(2,:,t) = mLogVols(1,:)';
		mStates(3,:,t) = mLogVols(2,:)';
		mWeights(:,t) = vWeights';
end

result.ESS = vESS;
result.states = mStates;
result.weights = mWeights;
result.filter = mFilteredEstimates;
result.predictor = mPredictedEstimates;
result.loglike = dLoglike;

return 

function [mLogVols, mFilter] = initializeParticles(iN)

phi = 0.99;
sigma2 = 0.05;
mu = -2.5;
mLogVols = zeros(2,iN);
mLogVols(:,:) = mu + sqrt(sigma2/(1-phi^2))*randn(2,iN);  % there is technically no stationary distribution
mFilter = zeros(2,iN);
mFilter(1,:) = 2.5 + sqrt(exp(mLogVols(1,:))).*randn(1,iN);
mFilter(2,:) = exp(mLogVols(2,:))*2;
return

function mLogVols = drawVolatilities(vParms, mLogVols, iN)
sigma2_meas = vParms(1,1);
sigma2_trans = vParms(2,1);
mLogVols(1,:) = mLogVols(1,:) + sqrt(sigma2_meas)*randn(1,iN); % simulate from the transition densities
mLogVols(2,:) = mLogVols(2,:) + sqrt(sigma2_trans)*randn(1,iN);
return

function mPredictions = KalmanPredictionStep(mLogVols, mFilter, iN)

	% Kalman filter matrices
	T = 1;
	c = 0;
	Q = exp(mLogVols(2,:));
	% One step ahead prediction
	mPredictions = zeros(2,iN);
	mPredictions(1,:) = T*mFilter(1,:) + c; 
    mPredictions(2,:) = T*mFilter(2,:)*T' + Q; 
return

function mFilter = KalmanUpdateStep(mLogVols, mPredictions, iN, vY)
	% Kalman filter matrices
	Z = 1;
    T = 1;
	H = exp(mLogVols(1,:));
    
    % prediction error
    v = vY - Z*mPredictions(1,:);
    % prediction error variance
	F = Z*mPredictions(2,:)*Z' + H;  
	% Kalman Gain
    K = T*mPredictions(2,:)*Z'.*(1 ./F);
    % Perform the Update Step
    mFilter = zeros(2,iN);
    mFilter(1,:) = mPredictions(1,:) + K.*v;
    mFilter(2,:) = mPredictions(2,:) - K.*Z.*mPredictions(2,:); 	
return

function vLogWeights = logIncrementalImportanceWeights(mLogVols, mPredictions, vY)

	% Kalman filter matrices
	Z = 1;
	H = exp(mLogVols(1,:));
	
	% prediction errors and prediction error variances
	v = vY - Z*mPredictions(1,:);
    F = Z*mPredictions(2,:)*Z' + H;  
   	vLogWeights = -0.5*log(2*pi)-0.5*log(abs(F)) -0.5*v.*(1 ./F).*v ;
return



function outIndex = systematicR_matlab(inIndex,wn)

if nargin < 2, error('Not enough input arguments.'); end

wn=wn';
[arb,N] = size(wn);  % N = Number of particles.

N_children=zeros(1,N);
label=zeros(1,N);
label=1:1:N;

s=1/N;
auxw=0;
auxl=0;
li=0;   % Label of the current point
% Initialisation
T=s*rand(1);
j=1;
Q=0;
i=0;

% Sampling before
u=rand(1,N);
while (T<1)
   if (Q>T)
      T=T+s;
      N_children(1,li)=N_children(1,li)+1;
   else
      % select i uniformly between j and N
      i=fix((N-j+1)*u(1,j))+j;
      % save the associate characteristic
      auxw=wn(1,i);
      li=label(1,i);
      % update the cfd
      Q=Q+auxw;
      % swap 
      wn(1,i)=wn(1,j);
      label(1,i)=label(1,j);
      %wn(1,j)=auxw;
      %label(1,j)=li;
      j=j+1;
   end
end

% COPY RESAMPLED TRAJECTORIES:  
% ============================
index=1;
for i=1:N
  if (N_children(1,i)>0)
    for j=index:index+N_children(1,i)-1
      outIndex(j) = inIndex(i);
    end;
  end;   
  index= index+N_children(1,i);   
end
	
