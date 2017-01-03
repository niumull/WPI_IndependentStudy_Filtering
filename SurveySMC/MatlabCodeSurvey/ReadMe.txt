This file contains Matlab code that replicates some of the results in the paper ``A survey of sequential Monte Carlo methods for economics and finance,'' Econometric Reviews, Vol. 31, No. 3, pp. 245-296, 2012.

The code is free for academic and educational use. 

The code contains MEX written in C that implement different types of resampling algorithms: 
1. systematic resampling (Carpenter Clifford Fearnhead 1999 IEEE)
2. stratified resampling (Kitagawa 1996 JCGS)
3. multinomial resampling (original bootstrap filter)
4. residual resampling (Liu Chen 1998 JASA)
 
These algorithms can be used in future work to speed up a particle filter. The .mexw64 files were compiled for a 64 bit computer and will not run on 32 bit machines. 
Instead, you can use the resampling algorithms that have been coded in Matlab.

If you use the code (e.g. the resampling algorithms in the MEX files) in future work, I would appreciate you citing the paper.

The actual graphs in the paper were produced using the Ox code. 