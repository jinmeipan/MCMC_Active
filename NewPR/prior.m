%snowpit structure
classdef prior

properties
    
    %prior types
    prior_type
    
    %inputted snowpits, sps used for prior calc
    spsr
    
    %inpuuted SWE, available using generic prior
    month
    site
    SWEr
    snow_class
    
    %SWE ratio, fixed
    SWE_std_ratio
    
    %basic
    swe_mean
    swe_std
    sd_mean       %all sd prior is calculate from swe and average_density
    sd_std

    
    density_mean  %1-6 layers
    density_std   %1-6 layers
    T_mean        %1-6 layers, used 274-T
    T_std         %1-6 layers, used 274-T
    dmax_mean     %1-6 layers
    dmax_std      %1-6 layers
    pex_mean      %1-6 layers
    pex_std       %1-6 layers
    dz_mean       %1-6 layers
    dz_std        %1-6 layers
    
    soilT_mean       %used 279-T
    soilT_std        %used 279-T
    roughness_mean
    roughness_std
    mv_soil_mean
    mv_soil_std
    
    
    %lognormal parameters, calculate directly from mean and std.
    %note that the values were checked to be not far from the values fitted
    %using lognormal distribution
    swe_mu
    swe_sigma
    sd_mu
    sd_sigma
    
    density_mu    %1-6 layers
    density_sigma %1-6 layers
    T_mu          %1-6 layers
    T_sigma       %1-6 layers
    dmax_mu       %1-6 layers
    dmax_sigma    %1-6 layers
    pex_mu        %1-6 layers
    pex_sigma     %1-6 layers
    dz_mu         %1-6 layers
    dz_sigma      %1-6 layers

	soilT_mu
    soilT_sigma
    roughness_mu
    roughness_sigma
    mv_soil_mu
    mv_soil_sigma
    
    %note
    note
end %properties


methods
    function prior=process(prior)
        
        if(strcmp(prior.prior_type,'local')==1)
            %first step, relayering the snowpits from 1 to 6 layers
            prior=relayering(prior);

            %second step, calculate mean and std
            prior=calc_local(prior);
        else            
            %zero step, read SWE
            prior=read_swe(prior);
            
            %first step, read snowtypes
            prior=read_class(prior);
        end
        
        prior = default_soil_mean(prior);
        prior = process_dz(prior);
        
        prior=calc_lognparam(prior);
    end
    
    
    
    
    function prior=reprocess(prior)
        %reprocess SWE, SD and dz prior given revised density prior, or SWE_std_ratio
        
        %reprocess SWE
        if(isnan(prior.SWE_std_ratio)==0)
            prior.swe_std=prior.swe_mean*prior.SWE_std_ratio;
        end
        
        %reprocess SD and dz
        prior=process_dz(prior);
        
        
        %recalc lognmal parameters
        prior=calc_lognparam(prior);
    end
    
    
    
    

    

    
    
    %% for local prior
    function prior=relayering(prior)
        for ip=1:length(prior.spsr)
            
            snowpit=prior.spsr(ip);
            
            for ilyr=1:6
                %reread Y, because Sodankyla snowpit has adde ice lense
                %information
                Zeros=zeros(snowpit.nlayer,1);
                %note that, I added one more field as dmax than what is
                %required by MEMLS
                Y_in = [[1:1:snowpit.nlayer]', snowpit.T, Zeros, snowpit.density, snowpit.dz*100, Zeros, snowpit.pex, snowpit.dmax];
                
                Y_out=resample2(Y_in,ilyr);
                
                snowpit.Y{ilyr+1}=Y_out;
                
                
%                 disp('Y_in');
%                 disp(num2str(Y_in));
%                 disp('Y_out');
%                 disp(num2str(Y_out));
                
            end
            
            prior.spsr(ip)=snowpit;
        end
    end
   
    
    
    function prior=calc_local(prior)
        prior=extract_prop(prior,'SWE');
        
        prior=extract_prop(prior,'soilT');
        %prior=extract_prop(prior,'roughness_fitted');  %because there is no roughness measurements
        prior=extract_prop(prior,'mv_soil');
        
        prior=extract_prop2(prior,'density');
        prior=extract_prop2(prior,'T');
        prior=extract_prop2(prior,'dmax');
        prior=extract_prop2(prior,'pex');
    end
    
    
    
    function prior=extract_prop(prior,value)
        
        np = length(prior.spsr);
        data=nan(np,1);
        
        for i=1:np
            strr = ['data(i) = prior.spsr(i).',value,';'];
            eval(strr);
        end
        
        idx=find(isnan(data));
        data(idx)=[];
        
        %
        if(strcmp(value,'SWE')==1)
            value='swe';
        end
        
        means=mean(data);
        stds=std(data);
        
        eval(['prior.',value,'_mean=means;'])
        eval(['prior.',value,'_std=stds;'])
    end
    
    
    
    function prior=extract_prop2(prior,value)
        
        %Y_in = [[1:1:snowpit.nlayer]', snowpit.T, Zeros, snowpit.density, snowpit.dz*100, Zeros, snowpit.pex, snowpit.dmax];   
        switch value
            case 'density'
                iss=4;
            case 'T'
                iss=2;
            case 'dmax'
                iss=8;
            case 'pex'
                iss=7;
        end
        
        np = length(prior.spsr);
        
        for ilyr=1:6
            eval(['prior.',value,'_mean{ilyr,1}=nan(ilyr,1);'])
            eval(['prior.',value,'_mean{ilyr,1}=nan(ilyr,1);'])
            
            for i=1:ilyr
                data = nan(np,1);
                for ip=1:np
                    data(ip)= prior.spsr(ip).Y{ilyr+1}(i,iss);
                end
                
                %summarize
                idx=find(isnan(data));
                data(idx)=[];
                
                means=mean(data);
                stds=std(data);
        
                eval(['prior.',value,'_mean{ilyr,1}(i,1)=means;'])
                eval(['prior.',value,'_std{ilyr,1}(i,1)=stds;'])
            end
        end
        
    end    
    
    
    
    %% for generic prior
    function prior=read_swe(prior)
        %set SWEr
        %extract SWE, accoring to month-site
        switch prior.site
            case 'sdkl'
              SWEr=[1	65.57	52.13293754
                    2	89.02	58.48861444
                    3	115.45	58.9826859
                    4	119.75	66.13780491
                    5	51.95	52.3943994
                    6	1.00	2.46178372
                    7	0.00	0.001221018
                    8	0.01	0.012110279
                    9	0.28	0.513010255
                    10	6.94	6.546879798
                    11	23.01	17.93403656
                    12	41.17	34.91444238];
            case 'churchill'
              SWEr=[1	77.20	25.47
                    2	95.79	24.41
                    3	108.37	25.69
                    4	104.59	35.78
                    5	45.7	40.9
                    6	2.3	8.2
                    7	0.0	0.0
                    8	0.0	0.0
                    9	0.4	0.5
                    10	4.7	4.2
                    11	28.8	14.3
                    12	60.0	24.0];
            case 'clpx'
              SWEr=[1	66.355	nan
                    2	90.522	38.8
                    3	109.773	38.8 %copied from Feb
                    4	105.267	nan
                    5	33.997	nan
                    6	0	nan
                    7	0	nan
                    8	0	nan
                    9	0	nan
                    10	0.819	nan
                    11	10.24	nan
                    12	35.226	nan];
        end
        
        prior.swe_mean=SWEr(prior.month,2);
        prior.swe_std=SWEr(prior.month,3);
    end
    
    
    
    function prior=read_class(prior)
        %set according to snow_type
        %determine density, 274-T, dmax and pex
        switch prior.snow_class
            case 'taiga'
                T_means=-11;
                T_stds=11;
                density_means= 217;
                density_stds=  56;
            case 'tundra'
                T_means=-18;
                T_stds=12;
                density_means= 284;
                density_stds=  75;
            case 'alpine'
                T_means=-8;
                T_stds=6;
                density_means= 335;
                density_stds=  86;
            case 'maritime'
                T_means=-5;
                T_stds=5;
                density_means= 343;
                density_stds= 101;
            case 'prairie'
                T_means=nan; %unknown
                T_stds=nan;
                density_means=312;
                density_stds = 85;
        end
        T_means=T_means+273.15;
        
        
        for ilyr=1:6
            prior.T_mean{ilyr}=ones(ilyr,1) * (T_means);
            prior.T_std{ilyr}=ones(ilyr,1) * T_stds;
            
            prior.density_mean{ilyr}=ones(ilyr,1) * density_means;
            prior.density_std{ilyr}=ones(ilyr,1) * density_stds;
        
            %the other parameters, fixed dmax and pex
            prior.dmax_mean{ilyr}=ones(ilyr,1)* 1;
            prior.dmax_std{ilyr}=ones(ilyr,1)* 1;
            prior.pex_mean{ilyr}=ones(ilyr,1)* 0.18;
            prior.pex_std{ilyr}=ones(ilyr,1)* 0.18;
        end
    end
            
    
    
    
    
    %% secondary functions
    function prior = default_soil_mean(prior)
        if(isempty(prior.soilT_mean)==1)
            %as the bottom-most soil temperature
            prior.soilT_mean = prior.T_mean{6}(1);
            prior.soilT_std = prior.T_std{6}(1);
        else if (isnan(prior.mv_soil_mean)==1)
                prior.soilT_mean = prior.T_mean{6}(1);
                prior.soilT_std = prior.T_std{6}(1);
            end
        end
        
        if(isempty(prior.roughness_mean)==1)
            prior.roughness_mean=1/100; %m
            prior.roughness_std=1/100;
        else if (isnan(prior.roughness_mean)==1)
                 prior.roughness_mean=1/100; %m
                 prior.roughness_std=1/100;
            end
        end
        
        if(isempty(prior.mv_soil_mean)==1)
            prior.mv_soil_mean=0.08; %m3/m3
            prior.mv_soil_std=0.08;
        else if (isnan(prior.mv_soil_mean)==1)
                prior.mv_soil_mean=0.08; %m3/m3
                prior.mv_soil_std=0.08;
            end
        end
    end
    
    
    function prior=process_dz(prior)
        %reprocess SD
        SWE_mean=prior.swe_mean;         %mm
        SWE_std=prior.swe_std;           
        density_mean=prior.density_mean{1}; %kg/m^3
        density_std=prior.density_std{1};   %kg/m^3
        
        %revised, add correlaction
        sd_mean=SWE_mean/density_mean;
        
        corr=0.3;
        sd_cov=SWE_mean^2./density_mean^2 * (SWE_std^2/SWE_mean^2 + density_std^2/...
            density_mean^2 - 2*corr*SWE_std*density_std / (SWE_mean*density_mean) );
        prior.sd_mean=sd_mean;
        prior.sd_std=sd_cov^0.5;
        
        %reprocess dz
        for i=1:6
            mean0=sd_mean/i;
            std0=sqrt(sd_cov/i);
            prior.dz_mean{i,1}=ones(i,1)*mean0;
            prior.dz_std{i,1}=ones(i,1)*std0;
        end
    end
    
    
    
    function prior=calc_lognparam(prior)
        
        m=prior.swe_mean;
        v=prior.swe_std.^2;
        prior.swe_mu = log((m^2)/sqrt(v+m^2));
        prior.swe_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.sd_mean;
        v=prior.sd_std.^2;
        prior.sd_mu = log((m^2)/sqrt(v+m^2));
        prior.sd_sigma = sqrt(log(v/(m^2)+1));
        
        
        %for each layers
        for i=1:6
            m=prior.density_mean{i};
            v=prior.density_std{i}.^2;
            prior.density_mu{i,1} = log((m.^2)./sqrt(v+m.^2));
            prior.density_sigma{i,1} = sqrt(log(v./(m.^2)+1));

            %m=prior.T_mean{i};
            m=274-prior.T_mean{i};    
            v=prior.T_std{i}.^2;
            prior.T_mu{i,1} = log((m.^2)./sqrt(v+m.^2));
            prior.T_sigma{i,1} = sqrt(log(v./(m.^2)+1));

            m=prior.dmax_mean{i};
            v=prior.dmax_std{i}.^2;
            prior.dmax_mu{i,1} = log((m.^2)/sqrt(v+m.^2));
            prior.dmax_sigma{i,1} = sqrt(log(v./(m.^2)+1));

            m=prior.pex_mean{i};
            v=prior.pex_std{i}.^2;
            prior.pex_mu{i,1} = log((m.^2)./sqrt(v+m.^2));
            prior.pex_sigma{i,1} = sqrt(log(v./(m.^2)+1));

            m=prior.dz_mean{i};
            v=prior.dz_std{i}.^2;
            prior.dz_mu{i,1} = log((m.^2)./sqrt(v+m.^2));
            prior.dz_sigma{i,1} = sqrt(log(v./(m.^2)+1));
        end
        
        %for soil parameters
        %m=prior.soilT_mean;
        m=279-prior.soilT_mean;
        v=prior.soilT_std.^2;
        prior.soilT_mu = log((m^2)/sqrt(v+m^2));
        prior.soilT_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.roughness_mean;
        v=prior.roughness_std.^2;
        prior.roughness_mu = log((m^2)/sqrt(v+m^2));
        prior.roughness_sigma = sqrt(log(v/(m^2)+1));
        
        m=prior.mv_soil_mean;
        v=prior.mv_soil_std.^2;
        prior.mv_soil_mu = log((m^2)/sqrt(v+m^2));
        prior.mv_soil_sigma = sqrt(log(v/(m^2)+1));
    end
   
end  %methods
end  %classdef
