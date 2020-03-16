# EarlyLife
Companion for "Turbulent coherent structures and early life below the Kolmogorov scale

System requirements
-------------------
The code was written in MATLAB version R2018a on a laptop running Mac OSX Mojave, and then run on the Harvard FAS Odyssey computing cluster to generate the results shown in the paper. However, the code uses very basic MATLAB functionalities and is therefore expected to run on most any version. No non-standard hardware is required; however, since the simulations are quite time-intensive, it would take many thousands of hours on a personal computer to achieve the number of realizations discussed in the manuscript. 


Installation guide
------------------
The code can run from the .m files included in this repository without any further changes to one's MATLAB installation. 



Demo and Instructions for Use -- main code, 'EarlyLife_SampleCode.m'
--------------------------------------------------------------------
The code can be called from the MATLAB command line with a small number of arguments. These are annotated in the included .m file. Expected run times can vary greatly, since a typical realization involves either extinction of a small inoculum of replicators (a fast process, which can be on the order of seconds to minutes) or the event we define as ``fixation'', meaning the population size increases 1000-fold, which can take several hours on a standard machine. 

Example call: 

tic
[outcome,tt]=EarlyLife_SampleCode(0.03,10,1,2)
toc

Note that if Replicase R1 is employed the simulations will necessarily be much slower, because the population size explodes much faster 


Demo and Instructions for Use -- code for Figure 4, 'WrightFisher_SampleCode.m'
--------------------------------------------------------------------
The code can be called from the MATLAB command line with a small number of arguments. These are annotated in the included .m file. This code should run extremely quickly relative to the main code.


Demo and Instructions for Use -- code for Figure 3, 'OneVort_SampleCode.m'
--------------------------------------------------------------------
The code can be called from the MATLAB command line with a small number of arguments. These are annotated in the included .m file. This code simply simulates the motion of two particles in the unsteady-vortex flow. 
