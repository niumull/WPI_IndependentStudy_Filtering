function vIndex = residualResampling(vWeights)

% Use Matlab's random number generator to generate the uniform r.v.'s
iN = size(vWeights,2);
% residual number of particles to sample
iN_residuals = iN - sum(floor(iN.*vWeights));
vU = rand(1,iN_residuals);

% Call the mex code to run in C
vIndex = residualResampling_C(vWeights, vU);

% Matlab implementation by Nando de Freitas and Arnaud for comparison
% vIndex = residualR_compare(1:iN, vWeights', vU);

% This C code uses an internal random number generator coded in C
% There is no need to pass the uniforms to it.
%vIndex = residualResampling_iC(vWeights);