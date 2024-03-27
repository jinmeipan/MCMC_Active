1. Open MCMC_in.m, revise paths, parameters and run
This code will generate texts for priors and observations for each snowpit, to be loaded by the MCMC fortran package.

2. Run fortran package in MCMC_V3_revised2_new2_normal_prior
The main function of this Fortran package is MetroMEMS.f90
This package will read inputs generated in 1 (FILENAME, RunParams.txt, tb_obs.txt, hyperpar.txt), and give three MCMC outputs (post_tb.out, post_theta.out, acceptance.txt)

3. Open MCMC_out.m, revise paths, and read the results from the Fortran package
This code will extract MCMC outputs from the Fortran folders and summarize the results of all snowpits.

Note:
1. Examples of the results are shown in the Active_Real_0.5_Mv&Sig_together folder.
2. There are several .mat files containing Matlab-defined classes. Add the path of the main folder and ./common_codes/functions before loading them.
