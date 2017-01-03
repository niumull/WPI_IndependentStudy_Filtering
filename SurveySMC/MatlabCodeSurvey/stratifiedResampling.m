function vIndex = stratifiedResampling(vWeights)

% Use Matlab's random number generator to generate the uniform r.v.'s
vU = rand(1,size(vWeights,2));

% Call the mex file
vIndex = stratifiedResampling_C(vWeights,vU);
