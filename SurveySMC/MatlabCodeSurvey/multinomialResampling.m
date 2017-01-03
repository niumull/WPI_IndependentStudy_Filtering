function vIndex = multinomialResampling(vWeights)

% Use Matlab's random number generator to generate the uniform r.v.'s
vU = rand(1,size(vWeights,2));

% Call the mex file to run the C code
vIndex = multinomialResampling_C(vWeights, vU);

% Matlab implementation by Nando de Freitas for comparison
% vIndex = multinomialR_compare(1:iN, vWeights', vU);

% This C code uses an internal random number generator coded in C
% vIndex = multinomialResampling_iC(vWeights);