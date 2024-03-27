% Prepare all MCMC Fortran code input files
% note, there are user-specific settings for priors and running options, in this main code and more in the functions!

clear all

%% set folders
package_folder='D:\Desktop\MCMC_Active-BASE-AM\';

%folder for code
folder_lib=package_folder;
%folder to output Fortran inputs
folder_out=[package_folder,'Active_Real_0.5_Mv&Sig_together\'];
%folder to retrieve prior information
folder_in=[package_folder,'NewPR\'];


%% start
addpath(folder_lib);
addpath([folder_lib,'common_codes'])
addpath([folder_lib,'common_codes\functions'])
addpath([folder_lib,'AtmosphericModel'])
cd(folder_lib);

%
tsp_file='sd_allsp.mat';  %matfile saved all the sodankyla snowpit measurements
prL_file='sd_L_other.mat'; %local prior from sodankyla monthly average
prG_file='sd_G.mat'; %generic prior from VIC
site='sd';
Opt_relative_dz=1;    %Jinmei,2017/06/30,relative thickness
Opt_Sigma_Error=0.5;  %Sigma Error, dB


%% produce output~
load([folder_in,tsp_file]); load([folder_in,prL_file]);  load([folder_in,prG_file]);
mkdir(folder_out);
cd(folder_out)


for ipits=1:69
    
    %% load pr
    sp=sp_tsp(ipits);                       
    pr=sd_G(sp.month);                    %to be revised!
    
    %% basic setting
    mcmc=MCMCRun4;
    mcmc.folder=[folder_out,site,num2str(ipits),'\'];  %subfolder for each pit~
    
    %backup snowpit information into mcmc
    mcmc.sp=sp;
    
    %set active backscattering frequency
    mcmc.active_freq=[10.2, 13.3, 16.7];  %to be revised!

    %set active observation angle
    mcmc.active_theta=[50];               %to be revised!
    
    %set active polarization
%     mcmc.active_pol={'vv','vh'};               %to be revised! {'vv','hh','vh'}; changed from vv to vh~
    mcmc.active_pol={'vv'};
    
    %set passive observations...
    mcmc.passive_freq=[];                 %to be revised!
    mcmc.passive_theta=[];                %to be revised!
    mcmc.passive_pol=[];                  %to be revised!
    
    %used to set other measurements, usually set as []
    mcmc.other_measurements=[];           %to be revised!
    
    %used to set snow class for your snowpit
    mcmc.SturmClass='taiga';              %to be revised!
    
    %total number of iterations
    mcmc.Niter=20000;

    %number of iteratios in burn-in period
    mcmc.Nburn=2000;

    mcmc.stdSigma=Opt_Sigma_Error;       % to be revised!!
    
    
    %% write basic setting
    mcmc=prep_runparam4(mcmc);
    
    %% write filenames
    mcmc=prep_filename4(mcmc);
    
    
    %% write microwave obserations
    mcmc=prep_tbobs4(mcmc,'real',ipits);
    
    
    %% write priors
    %std/mean for SWE prior
    pr.SWE_std_ratio=1;
    
    %prior for snow density, mean and std. Here used the values for taiga
    %snow in Sturm's class
    rho_mean=217; 
    rho_std=56;
    
    for i=1:6
        pr.density_mean{i}=repmat(rho_mean,i,1);
        pr.density_std{i}=repmat(rho_std,i,1);
        
        %prior for exponential correlation length
        pr.pex_mean{i}=repmat(0.18,i,1);
        pr.pex_std{i}=repmat(0.18/2,i,1);  %revised from 0.18 to 0.18/2 due to normal distr.
    end

    %temperature prior for snow
    for ilyr=1:6
        T_means=-10+273.15;
        T_stds=5;
        pr.T_mean{ilyr}=ones(ilyr,1) * (T_means);
        pr.T_std{ilyr}=ones(ilyr,1) * T_stds;
    end
    
    %temperature prior for soil
    pr.soilT_mean = -10+273.15;
    pr.soilT_std = 5;

    %soil moisture prior
    pr.mv_soil_mean = 0.08; %changed back to 8% volumetric content
    pr.mv_soil_std = 0.08/2;

    %soil roughness prior
    pr.roughness_mean = 1/100;  %1-cm roughness, to be revised, like 1-mm
    pr.roughness_std = (1/2)/100;
        
    %process and save prior 
    pr=reprocess(pr);
    prep_hyper4(mcmc,pr,Opt_relative_dz);

    save([mcmc.folder,'MCMC.mat'],'mcmc');
end
