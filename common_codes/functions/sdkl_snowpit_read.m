function [sps_date,sps_sd,sps_swe,sps_density,sps_dmax]=sdkl_snowpit_read(years)

%% snowpits
load('D:\Desktop\MCMC_Active-BASE-AM\NewPR\sd_allsp.mat')

sps_date=nan;
sps_sd=nan;
sps_swe=nan;
sps_density=nan;
sps_dmax=nan;

for i=1:length(sp_tsp)
    sp=sp_tsp(i);
    date0=datenum(sp.year,sp.month,sp.date,sp.time,0,0);
    if(date0>datenum(years(1),7,1) & date0<datenum(years(end),7,1))
        sps_date=[sps_date;date0];
        sps_sd=[sps_sd;sp.SD];
        sps_swe=[sps_swe;sp.SWE];
        sps_density=[sps_density;sp.avg_density];
        sps_dmax=[sps_dmax;sp.avg_dmax];
    end
end

sps_date(1)=[];
sps_sd(1)=[];
sps_swe(1)=[];
sps_density(1)=[];
sps_dmax(1)=[];

end
