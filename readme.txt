1.Open MCMC_in.m, revise parameters and run
This code will generate text information of prior and observations for each snowpit, to be loaded into MCMC fortran package.

2.Run fortran package in MCMC_V3_revised2_new2_normal_prior
The main function of the fortran package is MetroMEMS.f90
This package will generate three MCMC outputs: post_tb.out, post_theta.out, acceptance.txt

3.Open MCMC_out.m, read the result from fortran package
This code will extract MCMC outputs from Fortran folders,and summarizes results of all snowpits.
