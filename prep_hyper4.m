%write hyper.txt
%
% ---  Discussion on how to determin lognormal mean and std -----
% m and v are mean and varianc of the data X
%m = 1;
%v = 2;
% mu and sigma are mean and standard reviation of log(X)
%mu = log((m^2)/sqrt(v+m^2));
%sigma = sqrt(log(v/(m^2)+1));
%[exp(mu-sigma),exp(mu),exp(mu+sigma)]=
%0.202411359667921       0.577350269189626           1.6468113937884;
%[m-sqrt(v),m,m+sqrt(v)]=
%-0.414213562373095                         1          2.41421356237309;
%if the first run starts from exp(mu+sigma), than this method is
%acceptable
% ------------------------------------------------------------------


%thickness, temperature, grain size, density
function prep_hyper4(MCMCRun4,pr,Opt_relative_dz)


folder0=MCMCRun4.folder;


%% wrtie
fid_hyper=fopen([folder0,'hyperpar.txt'],'w');
Max_lyr=6;

for ilyr=1:Max_lyr
    
    fprintf(fid_hyper,'%2i\n',ilyr);
    
    %thickness
    fprintf(fid_hyper, 'Layer thickness (log-normal),m \n');
    
    if(Opt_relative_dz==0 | ilyr==1)
        for k=1:ilyr
            %fprintf(fid_hyper, '%12.8f\n',pr.dz_mu{ilyr}(k));
            fprintf(fid_hyper, '%12.8f\n',pr.dz_mean{ilyr}(k));
        end
        for k=1:ilyr
            %revised to normal distr. param
            fprintf(fid_hyper, '%12.8f\n',(pr.dz_std{ilyr}(k)/2).^2);
        end
    else
        %Jinmei, 2017/07/30, change to relative thickness as prior
        %it means the first layer will give absolute thickness; 
        %then, the following layers will give the relative thickness
        %compared to the first layer!
        m=1; v=0.2^2;
        mu0 = log((m^2)/sqrt(v+m^2));
        sigma0 = sqrt(log(v/(m^2)+1)); 
        
        %revised to normal distr. param
        fprintf(fid_hyper, '%12.8f\n',pr.dz_mean{ilyr}(k));
        for k=2:ilyr
            fprintf(fid_hyper, '%12.8f\n',m);
        end
        
        fprintf(fid_hyper, '%12.8f\n',(pr.dz_std{ilyr}(k)/2).^2);
        for k=2:ilyr
            fprintf(fid_hyper, '%12.8f\n',v);
        end
    end
    
    
    %grain size
    fprintf(fid_hyper, 'Grain size (log-normal),mm \n');
    
    for k=1:ilyr
        switch MCMCRun4.ModelOpt
            case 2
                fprintf(fid_hyper, '%12.8f\n',pr.dmax_mean{ilyr}(k));
            case 1
                fprintf(fid_hyper, '%12.8f\n',pr.pex_mean{ilyr}(k));
        end
    end
    for k=1:ilyr
        switch MCMCRun4.ModelOpt
            case 2
                fprintf(fid_hyper, '%12.8f\n',pr.dmax_std{ilyr}(k).^2);
            case 1
                fprintf(fid_hyper, '%12.8f\n',pr.pex_std{ilyr}(k).^2);
        end
    end
    
    
    %density
    fprintf(fid_hyper, 'Density (log-normal),kg/m^3 \n');
    
    for k=1:ilyr
        fprintf(fid_hyper, '%12.8f\n',pr.density_mean{ilyr}(k));
    end
    for k=1:ilyr
        fprintf(fid_hyper, '%12.8f\n',pr.density_std{ilyr}(k).^2);
    end
    
    
    %snow temperature
    fprintf(fid_hyper, '274.-Temperature (log-normal),K \n');
    
    for k=1:ilyr
        fprintf(fid_hyper, '%12.8f\n',274-pr.T_mean{ilyr}(k));  %snow T from bottom to surface
    end
    

    for k=1:ilyr
        fprintf(fid_hyper, '%12.8f\n',pr.T_std{ilyr}(k).^2);
    end
    
    
    %soil temperature
    fprintf(fid_hyper, 'Soil Temperature, 279-Tempereature (log-normal),K \n');
    fprintf(fid_hyper, '%12.8f\n',279-pr.soilT_mean);	 %ground T
    fprintf(fid_hyper, '%12.8f\n',pr.soilT_std.^2);
    
    
    %soil moisture and roughness
    fprintf(fid_hyper, 'Soil volumetric water content,FRAC \n');
%     fprintf(fid_hyper, '%12.8f\n',pr.mv_soil_mu);       %soil moisture
%     fprintf(fid_hyper, '%12.8f\n',pr.mv_soil_sigma^2);
    fprintf(fid_hyper, '%12.8f\n',pr.mv_soil_mean);       %soil moisture
    fprintf(fid_hyper, '%12.8f\n',pr.mv_soil_std^2);
    
    %soil roughness
	fprintf(fid_hyper, 'Soil surface rougness,m \n');
%     fprintf(fid_hyper, '%12.8f\n',pr.roughness_mu);       %soil roughness
%     fprintf(fid_hyper, '%12.8f\n',pr.roughness_sigma.^2);
    fprintf(fid_hyper, '%12.8f\n',pr.roughness_mean);
    fprintf(fid_hyper, '%12.8f\n',pr.roughness_std.^2);
    
    %model paramter, M
    fprintf(fid_hyper, 'Model Parameter with normal distribution, P_M \n');
    fprintf(fid_hyper, '%12.8f\n',0.1);
    fprintf(fid_hyper, '%12.8f\n',0.05^2);
%     m=0.1; 
%     %v=0.1;  %Jinmei, revised on 2017/6/30, correct to variance
%     v=0.05^2;
%     mu0 = log((m^2)/sqrt(v+m^2));
%     sigma0 = sqrt(log(v/(m^2)+1));  
%     fprintf(fid_hyper, '%12.8f\n',mu0);
%     fprintf(fid_hyper, '%12.8f\n',sigma0.^2);
    
    
    %model parameter,Q
    fprintf(fid_hyper, 'Model Parameter with uniform distribution, P_Q \n');
    fprintf(fid_hyper, '%12.8f\n',0.1);
    fprintf(fid_hyper, '%12.8f\n',(0.01)^2);
%     m=0.1; 
%     %v=0.02;  %Jinmei, revised on 2017/6/30, correct to variance
%     v=0.02^2;
%     mu0 = log((m^2)/sqrt(v+m^2));
%     sigma0 = sqrt(log(v/(m^2)+1));  
%     fprintf(fid_hyper, '%12.8f\n',mu0);
%     fprintf(fid_hyper, '%12.8f\n',sigma0.^2);
    
    %model paramter, SR
    fprintf(fid_hyper, 'Model Parameter with uniform distribution, P_SR \n');
    fprintf(fid_hyper, '%12.8f\n',1.0);
    fprintf(fid_hyper, '%12.8f\n',0.25^2);
%     m=1; 
%     %v=0.25;  %Jinmei, revised on 2017/6/30, correct to variance
%     v=0.25^2;
%     mu0 = log((m^2)/sqrt(v+m^2));
%     sigma0 = sqrt(log(v/(m^2)+1));  
%     fprintf(fid_hyper, '%12.8f\n',mu0);
%     fprintf(fid_hyper, '%12.8f\n',sigma0.^2);
end


fclose(fid_hyper);
end
