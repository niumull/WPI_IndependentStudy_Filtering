function result = bootstrapFilterStochasticVolatility(vY, vParms, iN)

[iDim,iT] = size(vY);

% ==============================
% Storage Matrices
% ==============================
vESS = zeros(1,iT);
vLogVolsFilter = zeros(1,iT);
vVolsFilter = zeros(1,iT);
vLogVolsPredictor = zeros(1,iT);

% ===============================
% Initialize the particles
% ===============================
vLogVols = initializeParticles(iN, vParms);
dLoglike = 0;

	
for t=1:iT
	
    % =================================================== 
	% 	Draw particles
	% ===================================================
	vLogVols = importanceDensity(vParms, vLogVols, iN);
	
    % =================================================== 
	% Estimate the state variable (one step ahead predicted)
    % ===================================================
	vLogVolsPredictor(t) = sum(vLogVols);
	
	% =================================================== 
	% 	Calculate importance weights
	% ===================================================
	vLogWeights = logIncrementalImportanceWeights(vY(t), vLogVols);
    
    % =================================================== 
	%  Stabilize the importance weight
	% ==================================================
    dMaxWeight = max(vLogWeights);
	vWeights = exp(vLogWeights - dMaxWeight);
    
    % =================================================== 
	% 	Calculate the log-likelihood
	% ===================================================
	dLoglike = dLoglike + log(mean(vWeights)) + dMaxWeight;
		 		
	% =================================================== 
	% Normalize the weights
	% ===================================================
	vWeights = vWeights/sum(vWeights);
	
    % =================================================== 
	% Estimate the state variable (filtered)
    % ===================================================
	vLogVolsFilter(t) = sum(vWeights.*vLogVols);
    vVolsFilter(t) = sum(vWeights.*exp(vLogVols/2));
		
	% =================================================== 
	% Compute the ESS
	% ===================================================
	vESS(t) = 1/sum(vWeights.^2);
		
	% =================================================== 
	% Resample and reset the weights
	% ===================================================
    % NOTE: If you get an error here it is likely that the MEX file which contains the resampling step is not working
        % correctly. You can always use the first resampling algorithm
        % which is coded in Matlab
        
	vIndex = systematicR_matlab(1:iN,vWeights');
    %vIndex = residualResampling(vWeights);
	%vIndex = systematicResampling(vWeights);
    %vIndex = multinomialResampling(vWeights);
    %vIndex = stratifiedResampling(vWeights);
    vLogVols = vLogVols(vIndex);
end

result.ESS = vESS;
result.loglike = dLoglike;
result.logFilter = vLogVolsFilter;
result.filter = vVolsFilter;
result.predictor = vLogVolsPredictor;

return 

function result = initializeParticles(iN, vParms)

phi = vParms(1);
sigma2 = vParms(2);
mu = vParms(3);

result = mu + sqrt(sigma2/(1-phi^2))*randn(1,iN);% stationary distribution of the transition equation

return

function result = importanceDensity(vParms, vLogVols, iN)

phi = vParms(1);
sigma2 = vParms(2);
mu = vParms(3);

result = mu + phi*(vLogVols - mu) + sqrt(sigma2)*randn(1,iN);
return;

function result = logIncrementalImportanceWeights(dY, vLogVols)
result = -0.5*log(2*pi)-0.5*vLogVols -0.5*(1 ./exp(vLogVols))*(dY^2);
return;

function outIndex = systematicR_matlab(inIndex,wn);

if nargin < 2, error('Not enough input arguments.'); end

wn=wn';
[arb,N] = size(wn);  % N = Number of particles.

% SYSTEMATIC RESAMPLING:
% ====================

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
	
