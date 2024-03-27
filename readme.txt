1.Open MCMC_in.m, revise paths, parameters and run
This code will generate text information of prior and observations for each snowpit, to be loaded into MCMC fortran package.

2.Run fortran package in MCMC_V3_revised2_new2_normal_prior
The main function of the fortran package is MetroMEMS.f90
This package will generate three MCMC outputs: post_tb.out, post_theta.out, acceptance.txt

3.Open MCMC_out.m, revise paths, and read the result from the fortran package
This code will extract MCMC outputs from Fortran folders, and summarizes results of all snowpits.

Note:
Examples of the results are shown in the Active_Real_0.5_Mv&Sig_together folder.
