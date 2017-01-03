function vIndex = systematicResampling(vWeights)

% Use Matlab's random number generator to generate the uniform r.v.'s
vU = rand(1,size(vWeights,2)+1);
vIndex = systematicResampling_C(vWeights, vU);

% This C code uses an internal random number generator coded in C
%vIndex = systematicResampling_iC(vWeights);