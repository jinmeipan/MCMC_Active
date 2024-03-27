% Read MCMC results

folder='D:\Desktop\MCMC_Active-BASE-AM\';

addpath(folder);
folder_in=[folder,'Active_Real_0.5_Mv&Sig_together\']

%run settings
isave=0;  %whether you want to overwrite the mcmc-result.mat in Fortran folders
ShowDetails=1; %whether you want to plot details of results for each snowpit

%id of snowpits to be checked
icheck=1:69;
N=max(icheck);


%define result files
mcmc_sd=mcmc_results;  %note:sd means sodankyla here

dz_res=nan(N,2);
dz_obsr=nan(N,2);
pex_res=nan(N,2);
pex_obsr=nan(N,2);
dz_pr=nan(N,2);
pex_pr=nan(N,2);


%% starting, read MCMC results
for i=1:N
    
    if(sum(i==icheck)==0)  %check the snowpits in icheck
        continue
    end
    
     if(1)   %process mcmc results
        mcmc=MCMCRun4;
        filename=[folder_in,'sd',num2str(i),'\MCMC.mat'];
        load(filename);

        %read MCMC results
        mcmc.folder=[folder_in,'sd',num2str(i),'\'];
        mcmc.Run_location=0;  %skip two 4 bytes, use 1; otherwise, use 0
        if(0)
            %let probability calculation to chooese layers
            mcmc.Lyrplan=1;
        else
            %fixed choice of layers, with mandantory layer choice as mcmc.Nlyr_choose
            mcmc.Lyrplan=2;mcmc.Nlyr_choose=2;
        end

        mcmc.SturmClass='taiga';
        mcmc=mcmc.Main;

        symbol=['pit',num2str(i)]; %set title of the plots
        close all;mcmc.PlotProfileCompare(1,symbol,33);

        if(ShowDetails==0)
            mcmc.PlotProfileCompare(1,symbol,33);
            mcmc.PlotObsChains(1,symbol);
            mcmc.PlotParamResult('dz',1,symbol);
            mcmc.PlotParamResult('D',1,symbol);
            mcmc.PlotParamResult('rho',1,symbol);
            mcmc.PlotParamResult('Tsnow',1,symbol);
            mcmc.PlotParamResult('Tsoil',1,symbol);
            mcmc.PlotParamResult('mvs',1,symbol);
            mcmc.PlotParamResult('sig',1,symbol);
        %     mcmc.PlotParamResult('sd',1,symbol);
        %     mcmc.PlotParamResult('swe',1,symbol);
        %     mcmc.PlotParamResult('P_M',1,symbol);
            mcmc.PlotParamResult('P_Q',1,symbol);
        %     mcmc.PlotParamResult('P_SR',1,symbol);
        end

        %save results
        if(isave==1)
            save([folder_in,'sd',num2str(i),'/MCMC-res.mat'],'mcmc');
        end

     else  %directly read mcmc results
         load([folder_in,'sd',num2str(i),'/MCMC-res.mat'],'mcmc');
     end
 
    disp(['Number of layers chosen at Pit ',num2str(i),' : ',num2str(mcmc.nHat)]);
    disp(['True SWE:',num2str(mcmc.sp.SWE),' - MCMC SWE:',num2str(mcmc.sweHat)]);
    
    
    mcmc_sd=getmcmc(mcmc_sd,mcmc,i,'sd');
    
    
    %this code valid only for 2-layer choice
    if(1)
        Zeros=repmat(0,mcmc.sp.nlayer,1);
        Y_in=[[1:1:mcmc.sp.nlayer]', mcmc.sp.T, Zeros, mcmc.sp.density, mcmc.sp.dz*100, Zeros, mcmc.sp.pex];
        Y_in=resample2(Y_in,2);
        
        dz_obsr(i,:)=Y_in(:,5);
        pex_obsr(i,:)=Y_in(:,7);
    else
        dz_obsr(i,:)=mcmc.sp.Y{1}(:,5)';
        pex_obsr(i,:)=mcmc.sp.Y{1}(:,7)';
    end
    
    dz_res(i,:)=mcmc.md_dz{2,1};
    pex_res(i,:)=mcmc.md_D{2,1};
    
    dz_pr(i,:)=mcmc.DzMean{1,2}';
    pex_pr(i,:)=mcmc.DMean{1,2}';
end



%% plot time-series snow depth and pex estimations
%this code valid only for 2-layer choice
figure;
set(gcf,'color','w','position',[0,0,1200,300]);
subplot(1,2,1);
title('Bottom Layer')
plot(mcmc_sd.date,dz_obsr(:,1)/100,'bo'); hold on;
plot(mcmc_sd.date,dz_res(:,1),'ro');
plot(mcmc_sd.date,dz_pr(:,1),'g+');
legend('Obsr','MCMC','Prior');
xlabel('PitNo.');ylabel('Layer thickness (m)')

subplot(1,2,2);
title('Surface Layer')
plot(mcmc_sd.date,dz_obsr(:,2)/100,'bx'); hold on;
plot(mcmc_sd.date,dz_res(:,2),'rx');
plot(mcmc_sd.date,dz_pr(:,2).*dz_pr(:,1),'g+');
legend('Obsr','MCMC','Prior');
xlabel('PitNo.');ylabel('Layer thickness (m)')

%
figure;
set(gcf,'color','w','position',[0,0,1200,300]);
subplot(1,2,1);
title('Bottom Layer')
plot(mcmc_sd.date,pex_obsr(:,1),'bo'); hold on;
plot(mcmc_sd.date,pex_res(:,1),'ro');
plot(mcmc_sd.date,pex_pr(:,1),'g+');
legend('Obsr','MCMC','Prior');
xlabel('PitNo.');ylabel('P_{ec} (mm)')

subplot(1,2,2);
title('Surface Layer')
plot(mcmc_sd.date,pex_obsr(:,2),'bx'); hold on;
plot(mcmc_sd.date,pex_res(:,2),'rx');
plot(mcmc_sd.date,pex_pr(:,2),'g+');
legend('Obsr','MCMC','Prior');
xlabel('PitNo.');ylabel('P_{ec} (mm)')



%%
mcmc_sd.opt_recalc_swe=0;
mcmc_sd.opt_delete_outlier=1;
mcmc_sd=mcmc_sd.detect_outlier;



%% plots
PlotNo=2;
mcmc_sd.plot_obsr(88);
mcmc_sd.plot_obsr_error(89);

mcmc_sd.plot_scatter('sd',100,PlotNo);
mcmc_sd.plot_scatter('swe',101,PlotNo);
mcmc_sd.plot_scatter('rhoavg',102,PlotNo);
mcmc_sd.plot_scatter('Davg',103,PlotNo);

mcmc_sd.plot_timeseries('sd',201,PlotNo);
mcmc_sd.plot_timeseries('swe',202,PlotNo);
mcmc_sd.plot_timeseries('rhoavg',203,PlotNo);
mcmc_sd.plot_timeseries('Davg',204,PlotNo);
mcmc_sd.plot_timeseries('tsnowavg',205,PlotNo);
mcmc_sd.plot_timeseries('tsoil',206,PlotNo);
mcmc_sd.plot_timeseries('mvs',207,PlotNo);
mcmc_sd.plot_timeseries('sig',208,PlotNo);
mcmc_sd.plot_timeseries('Qavg',209,PlotNo);


%% export errors
fname=[folder_in,'!mcmc_res.txt'];
fid=fopen(fname,'a');
mcmc_sd.error_stat('sd',fid);
mcmc_sd.error_stat('swe',fid);
mcmc_sd.error_stat('rhoavg',fid);
mcmc_sd.error_stat('Davg',fid);
fclose(fid);
