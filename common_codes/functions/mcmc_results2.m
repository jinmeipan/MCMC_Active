%this is a convenient class used to save MCMC snow retrieval results
%each variable is a vector of Npits, where Npits is the total number of snowpits
%time, iop or location of the snowpits are also available
%variables are named as true_ or mcmc_
%several functions are defined to plot the 1:1 and time-series comparison
%of variables
%functions are also defined to calcualte the mean and std of variables,
%with an outlier detection.

%used for active radar retrival only!! if not, revise getmcmc function - measured and simulated observations


classdef mcmc_results2
    
properties
    %number of snowpits
    Npits
    Nf
    Np
    Ntheta
    date
    site
    site_readme
    opt_delete_outlier
    opt_recalc_swe
    
    %prior swe
    pr_swe
    pr_sd
    
    %mcmc snowpit properties
    mcmc_sd
    mcmc_swe
    mcmc_swe2
    mcmc_dz
    mcmc_rhoavg
    mcmc_rho
    mcmc_D
    mcmc_Davg
    mcmc_tsnowavg
    mcmc_tsoil
    mcmc_mvs
    mcmc_sig
    mcmc_Qavg
    
    %true snowpit properties
    true_sd
    true_swe
    true_dz
    true_rhoavg
    true_rho
    true_D
    true_Davg
    true_Davg_fitted
    true_tsnowavg
    true_tsoil
    true_mvs
    true_sig
    true_Qavg
    
    
    %outliers
    mcmc_sd_std
    mcmc_swe_std
    mcmc_Davg_std
    mcmc_rhoavg_std
    mcmc_nHat
    mcmc_outliers
    
    %observations
    mcmc_obsr
    true_obsr
    mcmc_obsr_simu %to reconfirm the observation model in MCMC
 
end
    
methods
    function mcmc_results=getmcmc(mcmc_results,mcmc,i,site)
       
        mcmc_results.Npits=max(mcmc_results.Npits,i);

        %read into snow properties using fields in mcmcresults
        %basic information
        mcmc_results.site(i,1)=str2num(mcmc.sp.site(end));
        
        switch site
            case 'sd'
                mcmc_results.site_readme={'1=iop1','2=iop2','3=iop3','4=iop4'};
        end
        
        mcmc_results.date(i,1)=datenum(mcmc.sp.year,mcmc.sp.month,mcmc.sp.date);
        
        %read prior
        mcmc_results.pr_sd(i,1)=mcmc.sdPr*100;
        mcmc_results.pr_swe(i,1)=mcmc.swePr;
        
        %read snow/soil properties
        mcmc_results.mcmc_sd(i,1)=mcmc.sdHat*100;   %cm
        mcmc_results.mcmc_swe(i,1)=mcmc.sweHat;
        mcmc_results.mcmc_swe2(i,1)=mcmc.sdHat*mcmc.RhoMean{1};
        mcmc_results.mcmc_dz{i,1}=mcmc.md_dz{mcmc.nHat};
        mcmc_results.mcmc_rho{i,1}=mcmc.md_rho{mcmc.nHat};
        mcmc_results.mcmc_rhoavg(i,1)=mcmc.rhoavgHat;
        mcmc_results.mcmc_D{i,1}=mcmc.md_D{mcmc.nHat};
        mcmc_results.mcmc_Davg(i,1)=mcmc.DavgHat;
        mcmc_results.mcmc_tsnowavg(i,1)=mcmc.TsnowavgHat;
        mcmc_results.mcmc_tsoil(i,1)=mcmc.md_Tsoil{mcmc.nHat};
        mcmc_results.mcmc_mvs(i,1)=mcmc.md_mvs{mcmc.nHat}*100;
        mcmc_results.mcmc_sig(i,1)=mcmc.md_sig{mcmc.nHat}*1000;
        mcmc_results.mcmc_Qavg(i,1)=mcmc.md_P_Q{mcmc.nHat};
        
        
        mcmc_results.true_sd(i,1)=mcmc.sp.SD*100;   %cm
        if(1)
            mcmc_results.true_swe(i,1)=mcmc.sp.SWE;
        else
%           sp.Y{1}=[[1:1:sp.nlayer]', sp.T, Zeros, sp.density, sp.dz*100, Zeros, sp.pex];
            Y00=mcmc.sp.Y{1};
            mcmc_results.true_swe(i,1)=sum(Y00(:,4).*Y00(:,5)/100);
        end
        
        
        %% revised starting Mar26
        iop=unique(mcmc_results.site);
        idx=find(iop==0);
        iop(idx)=[];
        
        load(['/Users/jinmeipan/Downloads/MCMC_DataPrep/codes_Read_Snowpit/sodankyla/sp_sdkl_iop',num2str(iop),'_V0.mat'])
        eval(['sps=sp_sdkl_iop',num2str(iop),';']);
        
        date0=datenum(mcmc.sp.year,mcmc.sp.month,mcmc.sp.date,mcmc.sp.time,0,0);
        
        dates0=0;
        for ipits=1:length(sps)
            dates0=[dates0;datenum(sps(ipits).year,sps(ipits).month,sps(ipits).date,sps(ipits).time,0,0)];
        end
        dates0(1)=[];
        
        diff=abs(date0-dates0); 
        idx=find(diff==min(diff));
        
        if(abs(dates0(idx)-date0)<=1)
            sp_tsp0=sps(idx);
            mcmc_results.true_dz{i,1}=sp_tsp0.dz;
            mcmc_results.true_rho{i,1}=sp_tsp0.density;
            mcmc_results.true_rhoavg(i,1)=sp_tsp0.avg_density;
            mcmc_results.true_D{i,1}=sp_tsp0.pex;
            mcmc_results.true_Davg(i,1)=sp_tsp0.avg_pex;   %to be revised, whether use sp_processed or not~
            mcmc_results.true_tsnowavg(i,1)=sp_tsp0.avg_T-273.15;
        else
            mcmc_results.true_dz{i,1}=NaN;
            mcmc_results.true_rho{i,1}=NaN;
            mcmc_results.true_rhoavg(i,1)=NaN;
            mcmc_results.true_D{i,1}=NaN;
            mcmc_results.true_Davg(i,1)=NaN;
            mcmc_results.true_tsnowavg(i,1)=NaN;
        end
        
        mcmc_results.true_tsoil(i,1)=mcmc.sp.soilT(1)-273.15;      %revised
        mcmc_results.true_mvs(i,1)=mcmc.sp.mv_soil(1)*100;  %percent
        
        mcmc_results.true_sig(i,1)=0.1/100*1000; %mm
        mcmc_results.true_Qavg(i,1)=0.1; %mm
        
        %outlier detection
        mcmc_results.mcmc_sd_std(i,1)=mcmc.sdStdHat*100;
        mcmc_results.mcmc_swe_std(i,1)=mcmc.sweStdHat;
        mcmc_results.mcmc_Davg_std(i,1)=mcmc.DavgStdHat;
        mcmc_results.mcmc_rhoavg_std(i,1)=mcmc.rhoavgStdHat;
        mcmc_results.mcmc_nHat(i,1)=mcmc.nHat;
        
        
        %measured and simulated observations
        if(mcmc.active_np>0 & mcmc.passive_np==0)
            obsr_sigma=nan;
            mcmc_sigma=nan;
            
            dum = mcmc.ObsPost(mcmc.Nburn+1:end,:,mcmc.nHat,1); %revise from i to 1
            Yhat = mean(dum); %this gives a row
            nf=length(mcmc.active_freq);
            ntheta=length(mcmc.active_theta);
            npol=mcmc.active_np;
            
            %Yhat = reshape(Yhat,nf,mcmc.active_np);
            
            Yhat2=nan(npol,ntheta,nf);
            for izz=1:npol*ntheta*nf
                Yhat2(izz)=Yhat(izz);
            end
            
            
            for itheta=1:ntheta
                for ipp=1:mcmc.active_np
                    for iff=1:nf
                        obsr_sigma=[obsr_sigma,mcmc.ObsMCMCIn(ipp,iff,itheta)];
                        mcmc_sigma=[mcmc_sigma,Yhat2(ipp,itheta,iff)];
                    end
                end
            end
            
            mcmc_results.mcmc_obsr(i,:)=mcmc_sigma(2:end);
            mcmc_results.true_obsr(i,:)=obsr_sigma(2:end);
            
            %mcmc_results.mcmc_obsr_postsimu
            
            mcmc_results.Nf=nf;
            mcmc_results.Np=mcmc.active_np;
            mcmc_results.Ntheta=length(mcmc.active_theta);
        else
%             error('mcmc_results deals only active measurements. revise function getmcmcm for passive measurements.')
            %!!!!!!!!
            %to be revised!!          !!!!
            %!!!!!!
            obsr_sigma=nan;
            mcmc_sigma=nan;
            
            dum = mcmc.ObsPost(mcmc.Nburn+1:end,:,mcmc.nHat,1); %revise from i to 1
            Yhat = mean(dum); %this gives a row
            nf=length(mcmc.passive_freq);
            ntheta=length(mcmc.passive_theta);
            npol=mcmc.passive_np;
            
            %Yhat = reshape(Yhat,nf,mcmc.active_np);
            
            Yhat2=nan(npol,ntheta,nf);
            for izz=1:npol*ntheta*nf
                Yhat2(izz)=Yhat(izz);
            end
            
            
            for itheta=1:ntheta
                for ipp=1:mcmc.active_np
                    for iff=1:nf
                        obsr_sigma=[obsr_sigma,mcmc.ObsMCMCIn(ipp,iff,itheta)];
                        mcmc_sigma=[mcmc_sigma,Yhat2(ipp,itheta,iff)];
                    end
                end
            end
            
            mcmc_results.mcmc_obsr(i,:)=mcmc_sigma(2:end);
            mcmc_results.true_obsr(i,:)=obsr_sigma(2:end);
            
            %mcmc_results.mcmc_obsr_postsimu
            
            mcmc_results.Nf=nf;
            mcmc_results.Np=mcmc.passive_np;
            mcmc_results.Ntheta=length(mcmc.passive_theta);
        end
    end
    
    
    %plot scatter plots
    function plot_scatter(mcmc_results,property,figureNo,plotNo)
        if(strcmp(property,'swe')==1 & mcmc_results.opt_recalc_swe==1)
            mcmc_results.mcmc_swe=mcmc_results.mcmc_swe2; %use the prior density to recalcualte swe instead
        end
        
        iops=unique(mcmc_results.site);
        idx=find(iops==0);
        iops(idx)=[];
        colorstr={'b','g','r','m'};
        
        figure(figureNo);hold on;
        for i=1:length(iops)
            idx=find(mcmc_results.site==iops(i));
            
            if(plotNo==1)
            eval(['plot(mcmc_results.true_',property,'(idx),mcmc_results.mcmc_',property,...
                '(idx),''bo'',''MarkerSize'',10,''color'',colorstr{i});']);
            else if(plotNo==2)
                  eval(['plot(mcmc_results.true_',property,'(idx),mcmc_results.mcmc_',property,...
                '(idx),''bx'',''MarkerSize'',10,''color'',colorstr{i});']); 
                else
                    eval(['plot(mcmc_results.true_',property,'(idx),mcmc_results.mcmc_',property,...
                '(idx),''b*'',''MarkerSize'',10,''color'',colorstr{i});']);
                end
            end
        end
        %legend('iop1','iop2','iop3','iop4','location','northwest');
        title(['iop',num2str(iops)]);
           
        %plot prior if available
        if(strcmp(property,'swe')==1 | strcmp(property,'sd')==1)
            eval(['plot(mcmc_results.true_',property,',mcmc_results.pr_',property,...
                ',''k.'',''MarkerSize'',10);']);
        end
        

        switch property
            case 'sd'
                xlabel('Measured SD (cm)');ylabel('MCMC SD (cm)');plot([0,150],[0,150],'k--');
            case 'swe'
                xlabel('Measured SWE (mm)');ylabel('MCMC SWE (mm)');plot([0,300],[0,300],'k--');
            case 'rhoavg'
                xlabel('Measured bulk density (kg/m^3)');ylabel('MCMC bulk density (kg/m^3)');
                plot([100,400],[100,400],'k--');
            case 'Davg'
                xlabel('Measured bulk p_{ex} (mm)');ylabel('MCMC bulk p_{ex} (mm)');
                plot([0.15,0.35],[0.15,0.35],'k--');
            case 'tsnowavg'
                xlabel('Measured bulk Tsnow (^oC)');ylabel('MCMC bulk Tsnow (^oC)');
                plot([-15,0],[-15,0],'k--');
            case 'tsoil'
                xlabel('Measured bulk Tsoil (^oC)');ylabel('MCMC bulk Tsoil (^oC)');
                plot([-15,0],[-15,0],'k--');
            case 'mvs'
                xlabel('Measured soil moisture (vol%)');ylabel('MCMC soil moisture (vol%)');
                plot([0,0.2],[0,0.2],'k--');
            case 'sig'
                xlabel('Fitted soil roughness (mm)');ylabel('MCMC soil roughness (mm)');
                plot([0,6],[0,6],'k--');
            case 'Qavg'
                xlabel('Fitted Q');ylabel('MCMC Q');
                plot([0.6,1.4],[0.6,1.4],'k--');
            otherwise
                warning(['axis title not avaible for ',property]) 
        end
                
    end
    
      
    function plot_timeseries(mcmc_results,property,figureNo,plotNo)
        if(strcmp(property,'swe')==1 & mcmc_results.opt_recalc_swe==1)
            mcmc_results.mcmc_swe=mcmc_results.mcmc_swe2; %use the prior density to recalcualte swe instead
        end
        
        
        fig=figure(figureNo);set(gcf,'Position',[0,0,600,300]);hold on;
        %set(fig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
    
        eval(['plot(mcmc_results.date,mcmc_results.true_',property,',''bo'',''LineWidth'',1.0);'])
        if(plotNo==1)
            eval(['plot(mcmc_results.date,mcmc_results.mcmc_',property,',''kx-'',''LineWidth'',1.0);'])
        else if(plotNo==2)
                eval(['plot(mcmc_results.date,mcmc_results.mcmc_',property,',''k^'',''LineWidth'',1.0);'])
            else
                eval(['plot(mcmc_results.date,mcmc_results.mcmc_',property,',''k+'',''LineWidth'',1.0);'])
            end
        end
        legend('measured','mcmc','location','northwest')
        
        if(strcmp(property,'swe')==1 | strcmp(property,'sd')==1)
            eval(['plot(mcmc_results.date,mcmc_results.pr_',property,',''g+'');'])
            legend('measured','mcmc','prior','location','northwest')
        end
    
%         if(strcmp(property,'Davg')==1)
%             plot(mcmc_results.date,mcmc_results.true_Davg_fitted,'md');
%             legend('measured','mcmc','fitted','location','northwest');
%         end

        iops=unique(mcmc_results.site);
        idx=find(iops==0);
        iops(idx)=[];
        title(['IOP',num2str(iops)],'Fontsize',16)
        
        years=iops(1)+2008;
        xlims=[datenum(years,10,1),datenum(years+1,7,1)];    
        set(gca,'xlim',[datenum(2009,12,1),datenum(2013,5,1)]);
        
        datetick('x','mmm');
        set(gca,'xlim',xlims);
    
        datetick('x','mmm','keeplimits');
        set(gca,'Ygrid','on');
    
    
        switch property
            case 'sd'
                ylabel('SD (cm)');
            case 'swe'
                ylabel('SWE (mm)');
            case 'rhoavg'
                ylabel('Bulk density (kg/m^3)');
            case 'Davg'
                ylabel('Bulk p_{ex} (mm)');
            case 'tsnowavg'
                ylabel('Bulk T_{snow} (^oC)');
            case 'tsoil' 
                ylabel('T_{soil} (^oC)');
            case 'mvs'
                ylabel('Soil moisture (vol%)')
            case 'sig'
                ylabel('Soil roughness (mm)')
            case 'Qavg'
                ylabel('Q')
            otherwise
                warning(['axis title not avaible for ',property]) 
        end

    end
    
    
    function plot_obsr(mcmc_results,figureNo)
        %mcmc_results.mcmc_obsr(i,1)=obsr_sigma(2:end);
        %mcmc_results.true_obsr(i,1)=mcmc_sigma(2:end);
        colorstr={'k','r','g'};
        ylabelstr={'\sigma_0^{vv}','\sigma_{vh}'};
        
        gcf=figure(figureNo);hold on;
        set(gcf,'Position',[0,0,600,300]);
        for itheta=1:mcmc_results.Ntheta
            thetaNo=(mcmc_results.Np)*(mcmc_results.Nf)*(itheta-1);
            
            for ipp=1:mcmc_results.Np
                
                subplot(mcmc_results.Ntheta,mcmc_results.Np, ipp + (itheta-1)*mcmc_results.Np); hold on;
                for iff=1:mcmc_results.Nf
                    plot(mcmc_results.date,mcmc_results.true_obsr(:,thetaNo+iff+(ipp-1)*mcmc_results.Nf),...
                        'ko','color',colorstr{iff},'LineWidth',1.0);
                    plot(mcmc_results.date,mcmc_results.mcmc_obsr(:,thetaNo+iff+(ipp-1)*mcmc_results.Nf),...
                        'kx','color',colorstr{iff},'LineWidth',1.0);
                end
    %             datetick('x','mm/dd/yy');
                ylabel(ylabelstr{ipp})
                
                if(mcmc_results.Np==1 & mcmc_results.Nf==3)
                    legend('10.2VV-measured','10.2VV-mcmc',...
                        '13.3VV-measured','13.3VV-mcmc',...
                        '16.7VV-measured','16.7VV-mcmc')
                end
        %         set(gca,'xlim',[datenum(2009,12,1),datenum(2013,5,1)]);

                iops=unique(mcmc_results.site);
                idx=find(iops==0);
                iops(idx)=[];
                years=iops(1)+2008;
                xlims=[datenum(years,10,1),datenum(years+1,7,1)];    
                set(gca,'xlim',[datenum(2009,12,1),datenum(2013,5,1)]);

                datetick('x','mmm');
                set(gca,'xlim',xlims);

                datetick('x','mmm','keeplimits');
                set(gca,'Ygrid','on');
            end
        end
    end   
    
    
    function plot_obsr_error(mcmc_results,figureNo)
        %mcmc_results.mcmc_obsr(i,1)=obsr_sigma(2:end);
        %mcmc_results.true_obsr(i,1)=mcmc_sigma(2:end);
        colorstr={'kx','ro','gsq'};
        ylabelstr={'\sigma_0^{vv}','\sigma_{vh}'};
        
        gcf=figure(figureNo);hold on;
        set(gcf,'Position',[0,0,600,300]);
        
        for itheta=1:mcmc_results.Ntheta
            thetaNo=(mcmc_results.Np)*(mcmc_results.Nf)*(itheta-1);
            
            for ipp=1:mcmc_results.Np
                subplot(mcmc_results.Ntheta,mcmc_results.Np, ipp + (itheta-1)*mcmc_results.Np); hold on;
                for iff=1:mcmc_results.Nf
    %                 plot(mcmc_results.date,mcmc_results.true_obsr(:,iff+(ipp-1)*mcmc_results.Nf),...
    %                     'ko','color',colorstr{iff});
    %                 plot(mcmc_results.date,mcmc_results.mcmc_obsr(:,iff+(ipp-1)*mcmc_results.Nf),...
    %                     'kx','color',colorstr{iff});
                     error_sigma=mcmc_results.mcmc_obsr(:,thetaNo+iff+(ipp-1)*mcmc_results.Nf) - ...
                         mcmc_results.true_obsr(:,thetaNo+iff+(ipp-1)*mcmc_results.Nf);
                     plot(mcmc_results.date,error_sigma,colorstr{iff},'LineWidth',1.0);
                end
                datetick('x','mm/dd/yy');
                ylabel(ylabelstr{ipp})
                
                % set(gca,'xlim',[datenum(2009,12,1),datenum(2013,5,1)]);
                if(mcmc_results.Np==1 & mcmc_results.Nf==3)
                    legend('10.2VV-error',...
                        '13.3VV-error',...
                        '16.7VV-error')
                end

                iops=unique(mcmc_results.site);        
                idx=find(iops==0);
                iops(idx)=[];
                years=iops(1)+2008;
                xlims=[datenum(years,10,1),datenum(years+1,7,1)];    
                set(gca,'xlim',[datenum(2009,12,1),datenum(2013,5,1)]);

                datetick('x','mmm');
                set(gca,'xlim',xlims);

                datetick('x','mmm','keeplimits');
                set(gca,'Ygrid','on');
            end
        end
    end  
    
      
    
    function mcmc_results=detect_outlier(mcmc_results)
        %method 1, use std to determine outliers
        %method 2, use the error of sigma simulation to detect the outliers
        if(mcmc_results.opt_recalc_swe==1)
            a=mcmc_results.mcmc_sd_std; %mcmc_sd_std;(?)
        else
            a=mcmc_results.mcmc_swe_std; %mcmc_sd_std;(?)
        end
        
        alpha=0.05;
        [b,idx,outliers]=deleteoutliers(a,alpha);
        
        mcmc_results.mcmc_outliers=a*0;
        mcmc_results.mcmc_outliers(idx)=1;
        
        figure(999);hold on;
        plot(mcmc_results.date,mcmc_results.true_swe,'bo');
        plot(mcmc_results.date,mcmc_results.mcmc_swe,'kx');
        plot(mcmc_results.date,a,'g^','MarkerSize',12);
        aa_axis=axis;
        for i=1:length(idx)
            plot([mcmc_results.date(idx(i)),mcmc_results.date(idx(i))],...
                 [aa_axis(3),aa_axis(4)],'g--','LineWidth',1.0);
        end
        legend('measured swe','mcmc swe','mcmc swe std','outliers')
    end

    
    
    
    function mcmc_results=error_stat(mcmc_results,property,fid)
        
        if(strcmp(property,'swe')==1 & mcmc_results.opt_recalc_swe==1)
            mcmc_results.mcmc_swe=mcmc_results.mcmc_swe2; %use the prior density to recalcualte swe instead
        end
        
        
        %for different iops
        iops=unique(mcmc_results.site);
        idx=find(iops==0);
        iops(idx)=[];
        disp(['MCMC error for ',property, ': mb rmse']);
        fprintf(fid,['MCMC error for ',property, ': mb rmse \n']);
        
        for i_iops=1:length(iops)
            if(mcmc_results.opt_delete_outlier==1) %exclude outliers
                idx=find(mcmc_results.mcmc_outliers==0 & mcmc_results.site==i_iops);
                eval(['err0=mcmc_results.mcmc_',property,'(idx)-mcmc_results.true_',property,'(idx);']);
                
                idx2=find(isnan(err0)==0); err0=err0(idx2);
                
                mb=mean(err0);
                rmse=sqrt(sum(err0.^2)/length(err0));
                
                disp(['iop',num2str(i_iops),':', num2str([mb,rmse])]);
                fprintf(fid,['iop',num2str(i_iops),':', num2str([mb,rmse]),'\n']);
            end
        end

        
        if(strcmp(property,'swe')==1 | strcmp(property,'sd')==1)
            
            disp(['Prior error for ',property, ': mb rmse']);
            fprintf(fid,['Prior error for ',property, ': mb rmse \n']);
            
            for i_iops=1:length(iops)
                if(mcmc_results.opt_delete_outlier==1) %exclude outliers
                    idx=find(mcmc_results.mcmc_outliers==0 & mcmc_results.site==i_iops);
                    eval(['err0=mcmc_results.pr_',property,'(idx)-mcmc_results.true_',property,'(idx);']);
                    
                    idx2=find(isnan(err0)==0); err0=err0(idx2);
                    
                    mb=mean(err0);
                    rmse=sqrt(sum(err0.^2)/length(err0));
                    
                    disp(['iop',num2str(i_iops),':', num2str([mb,rmse])]);
                    fprintf(fid,['iop',num2str(i_iops),':', num2str([mb,rmse]),'\n']);
                end
            end
        end
        
        if(strcmp(property,'rhoavg')==1)
            disp(['Prior error for ',property, ': mb rmse']);
            fprintf(fid,['Prior error for ',property, ': mb rmse \n']);
            
            for i_iops=1:length(iops)
                if(mcmc_results.opt_delete_outlier==1) %exclude outliers
                    idx=find(mcmc_results.mcmc_outliers==0 & mcmc_results.site==i_iops);
                    eval(['err0=217-mcmc_results.true_',property,'(idx);']);
                    
                    idx2=find(isnan(err0)==0); err0=err0(idx2);
                    
                    mb=mean(err0);
                    rmse=sqrt(sum(err0.^2)/length(err0));
                    
                    disp(['iop',num2str(i_iops),':', num2str([mb,rmse])]);
                    fprintf(fid,['iop',num2str(i_iops),':', num2str([mb,rmse]),'\n']);
                end
            end
        end
        
        
        if(strcmp(property,'Davg')==1)
            disp(['Prior error for ',property, ': mb rmse']);
            fprintf(fid,['Prior error for ',property, ': mb rmse \n']);
            
            for i_iops=1:length(iops)
                if(mcmc_results.opt_delete_outlier==1) %exclude outliers
                    idx=find(mcmc_results.mcmc_outliers==0 & mcmc_results.site==i_iops);
                    eval(['err0=0.18*1.2-mcmc_results.true_',property,'(idx);']);
                    
                    idx2=find(isnan(err0)==0); err0=err0(idx2);
                    
                    mb=mean(err0);
                    rmse=sqrt(sum(err0.^2)/length(err0));
                    
                    disp(['iop',num2str(i_iops),':', num2str([mb,rmse])]);
                    fprintf(fid,['iop',num2str(i_iops),':', num2str([mb,rmse]),'\n']);
                end
            end
        end
        
        
    end
end

end