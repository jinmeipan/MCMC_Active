%Run params prepare
%input:
%-model: which observation model; e.g. model='hut'; where 'emp'=empirical memls; 'iba'='iba','hut'=hut model combination
%-site: where is the pits; e.g. site='sdkl'
%-sp: the snowpit information in type of Structure Snowpit
%-stdTb: observation error of TB (K); eg. stdTb=2.0

function MCMCRun4=prep_runparam4(MCMCRun4)


folder0=MCMCRun4.folder;
mkdir(folder0);

%% basic settings
MCMCRun4.NlyrMax=2;
MCMCRun4.Npits=1;
MCMCRun4.stdTb=2;
                       % !!!!!!!!     important changes here!!
%MCMCRun4.stdSigma=0.5; %changed from 0.1 to 0.5 dB (hypothesis, 0.1 dB may result in overfit)
MCMCRun4.ModelOpt=1;   %1 (MEMLS), 2 (HUT), 3 (DMRT-ML), 4 (DMRT-QMS)
MCMCRun4.ScatOpt=2;    %2 (MEMLS-IBA)
MCMCRun.UsePrior=1;

%whether m,q,SR will be estimated in MCMC (parameters used for the MEMLS
%backscattering simulation)
Estimate_M=0;
Estimate_Q=1; %Estimated Q
Estimate_SR=0;

%combine frequencies and thetas
freq=[MCMCRun4.active_freq,MCMCRun4.passive_freq];
theta=[MCMCRun4.active_theta,MCMCRun4.passive_theta];

freq=unique(freq); theta=unique(theta);
freq=sort(freq); theta=sort(theta);
MCMCRun4.MCMC_freq=freq;
MCMCRun4.MCMC_theta=theta;

MCMCRun4.active_np=length(MCMCRun4.active_pol);
MCMCRun4.passive_np=length(MCMCRun4.passive_pol);
MCMCRun4.other_np=length(MCMCRun4.other_measurements);
MCMCRun4.MCMC_pol=[repmat('p',MCMCRun4.passive_np),repmat('a',MCMCRun4.active_np),repmat('o',MCMCRun4.other_np)];

%calculate atm Tb, using -15 degC, 80% RH, and 1000 hPa (no strong
%difference using other inputs)
Tair=-15+273.15;
pressure=1000;
rh=80;
moist=RH2AH(rh,Tair);
month=MCMCRun4.sp(1).month;

for i=1:length(freq)
    for j=1:length(theta)
        [Tsky, tran] = skytemp2(freq(i),theta(j),'U',0.05,pressure,Tair,month,moist,[0,0,0,0,0,0,0]);
        skytb(i,j)=Tsky;
    end
end
MCMCRun4.skytb=skytb;


%% begin to write
fid_rp=fopen([folder0,'RunParams.txt'],'w');

fprintf(fid_rp,'Number of pits to run (Npits)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.Npits);

fprintf(fid_rp,'Number of iterations in the Markov Chain (Niter)\n');
fprintf(fid_rp,'%12d\n',MCMCRun4.Niter);

fprintf(fid_rp,'Number of burn-in iterations in the Markov Chain (Nburn)\n');
fprintf(fid_rp,'%12d\n',MCMCRun4.Nburn);

fprintf(fid_rp,'Number of snow layers to predict (NlyrMax)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.NlyrMax);

fprintf(fid_rp,'Number of observation frequencies (Nfreq)\n');
fprintf(fid_rp,'%5d\n',length(freq));

fprintf(fid_rp,'Number of observation angles (Nangle)\n');
fprintf(fid_rp,'%5d\n',length(theta));

fprintf(fid_rp,'Number of polarizations of passive Tb (Np_passive)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.passive_np);

fprintf(fid_rp,'Number of polarizations of active backscattering coefficient (Np_active)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.active_np);

fprintf(fid_rp,'Number of other observations (Np_other)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.other_np);

fprintf(fid_rp,'Number of snow parameters per layer (NsnowVars)\n');
fprintf(fid_rp,'%5d\n',4);

fprintf(fid_rp,'Number of soil parameters (NsoilVars)\n');
fprintf(fid_rp,'%5d\n',3);

fprintf(fid_rp,'Scattering Coefficient Option (ScatOpt): 1 (Empirical MEMLS, Hallikainen-HUT), 2 (MEMLS-IBA, combined HUT), or 3 (Roy-HUT)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.ScatOpt);

fprintf(fid_rp,'Observation Model Option (ModelOpt): 1 (MEMLS), 2 (HUT), 3 (DMRT-ML), 4 (DMRT-QMS), 5 (Bi-Continuous)\n');
fprintf(fid_rp,'%5d\n',MCMCRun4.ModelOpt);

fprintf(fid_rp,'Error standard deviation of Tb observations (StdTb)\n');
fprintf(fid_rp,'%10.4f\n',MCMCRun4.stdTb);

fprintf(fid_rp,'Error standard deviation of backscattering coefficient observations (StdSigma)\n');
fprintf(fid_rp,'%10.4f\n',MCMCRun4.stdSigma);
disp(['The \sigma observation error is: ',num2str(MCMCRun4.stdSigma)])

fprintf(fid_rp,'Error standard deviation of other observations (StdOther(Np_other))(Here StdOther(1)=Snow Depth in m)\n');
if(MCMCRun4.other_np>0)
    fprintf(fid_rp,'%10.4f\n',MCMCRun4.stdOther);
end

fprintf(fid_rp,'Observation frequencies (Freq(Nfreq))\n');
for i=1:length(freq)
	fprintf(fid_rp,'%10.4f\n',freq(i));
end

fprintf(fid_rp,'Observation angles in degree (Angle(Nangle))\n');
for i=1:length(theta)
	fprintf(fid_rp,'%10.4f\n',theta(i));
end

fprintf(fid_rp,'Polarizations for Passive Measurement(Pol_Passive(Np_passive)): 1 (vertical), 2 (horizontal)\n');
if(MCMCRun4.passive_np>0)
    for i=1:MCMCRun4.passive_np
        switch MCMCRun4.passive_pol{i}
            case 'v'
                fprintf(fid_rp,'%5d\n',1);
            case 'h'
                fprintf(fid_rp,'%5d\n',2);
        end
    end
end

  
fprintf(fid_rp,'Polarizations for Active Measurement(Pol_Active(Np_active)): 1 (vv), 2 (hh), 3 (vh), 4 (hv)\n');
if(MCMCRun4.active_np>0)
    for i=1:MCMCRun4.active_np
        switch MCMCRun4.active_pol{i}
            case 'vv'
                fprintf(fid_rp,'%5d\n',1);
            case 'hh'
                fprintf(fid_rp,'%5d\n',2);
            case 'vh'
                fprintf(fid_rp,'%5d\n',3);
            case 'hv'
                fprintf(fid_rp,'%5d\n',3);
        end
    end
end
   
    
fprintf(fid_rp,'Tb boundary condition above snow surface at vertical and horizontal polarizations (K) (Tsky(Nfreq,2))\n');
for i=1:length(freq)
    for j=1:length(theta)
        fprintf(fid_rp,'%10.4f\n',skytb(i,j));
    end
end

   
fprintf(fid_rp,'Use prior information, 0 or 1. (UsePrior)\n');
fprintf(fid_rp,'%5d\n',MCMCRun.UsePrior);


fprintf(fid_rp,'Estimate dZ, -1, 0, 1  (EstimateDz)\n');
fprintf(fid_rp,'%5d\n',1);

fprintf(fid_rp,'Estimate rho, -1, 0, 1  (EstimateRho)\n');
fprintf(fid_rp,'%5d\n',1);

fprintf(fid_rp,'Estimate D, -1, 0, 1  (EstimateD)\n');
fprintf(fid_rp,'%5d\n',1);

fprintf(fid_rp,'Estimate T, -1, 0, 1  (EstimateTsnowTsoil)\n');
fprintf(fid_rp,'%5d\n',1);          %% note, this will determine whether to iteratively estimate
                                    %both snow temperature and soil
                                    %temperature 

fprintf(fid_rp,'Estimate S, -1, 0, 1  (EstimateSoilRoughness)\n');
fprintf(fid_rp,'%5d\n',1);           %% note, this will determine whether to estimate soil roughness,
                                     %whereas the soil moisture will always
                                     %estimate

fprintf(fid_rp,'Estimate P_M, -1,0,1 (EstimateP_M)\n');
fprintf(fid_rp,'%5d\n',Estimate_M);

fprintf(fid_rp,'Estimate P_Q, -1,0,1 (EstimateP_Q)\n');
fprintf(fid_rp,'%5d\n',Estimate_Q);

fprintf(fid_rp,'Estimate P_SR, -1,0,1 (EstimateP_SR)\n');
fprintf(fid_rp,'%5d\n',Estimate_SR);

fprintf(fid_rp,'Initial Value dZ, -1, 0, 1  (InitialDz): 0 (normal value, mu), -1 (smaller value, mu-std), 1 (larger value, mu+std)\n');
fprintf(fid_rp,'%5d\n',0);

fprintf(fid_rp,'Initial Value rho, -1, 0, 1  (InitialRho)\n');
fprintf(fid_rp,'%5d\n',0);

fprintf(fid_rp,'Initial Value D, -1, 0, 1  (InitialD)\n');
fprintf(fid_rp,'%5d\n',0);

fprintf(fid_rp,'Initial Value T, -1, 0, 1  (InitialTsnow)\n');
fprintf(fid_rp,'%5d\n',0);

fprintf(fid_rp,'Initial Value S, -1, 0, 1  (InitialSoil)\n');
fprintf(fid_rp,'%5d\n',0);

fprintf(fid_rp,'Constrain rho (ConstrainRho): 0 (No constraint), 1 (surface<bottom, loose snow at surface), 2 (surface>bottom, dense snow at surface)\n');
fprintf(fid_rp,'%5d\n',1);

fprintf(fid_rp,'Constrain 274-T (ConstrainTsnow): 0 (No constraint), 1 (surface<bottom, warm snow at surface), 2 (surface>bottom, cold snow at surface)\n');
fprintf(fid_rp,'%5d\n',2);


fprintf(fid_rp,'Minimum & maximum limits for layer thickness [m]\n');
fprintf(fid_rp,'%10.4f\n',0.001);
fprintf(fid_rp,'%10.4f\n',10.0);

fprintf(fid_rp,'Minimum & maximum limits for density [kg/m3]\n');
fprintf(fid_rp,'%10.4f\n',50.0);
fprintf(fid_rp,'%10.4f\n',917.0);

fprintf(fid_rp,'Minimum & maximum limits for grain diameter [mm]\n');
fprintf(fid_rp,'%10.4f\n',0.001); 
fprintf(fid_rp,'%10.4f\n',5.0);

fprintf(fid_rp,'Minimum & maximum limits for snow temperature [degC]\n');
fprintf(fid_rp,'%10.4f\n',-30.0);
fprintf(fid_rp,'%10.4f\n',0.00);

fprintf(fid_rp,'Minimum & maximum limits for soil moisture [frac]\n');
fprintf(fid_rp,'%10.4f\n',0.0);
fprintf(fid_rp,'%10.4f\n',1.0);

fprintf(fid_rp,'Minimum & maximum limits for soil rms-height [m]\n');
fprintf(fid_rp,'%10.4f\n',0.0);
fprintf(fid_rp,'%10.4f\n',0.1);  %revised maximum of GndSig

fprintf(fid_rp,'Minimum & maximum limits for soil temperature [degC]\n');
fprintf(fid_rp,'%10.4f\n',-30.0);
fprintf(fid_rp,'%10.4f\n',5.0);

fprintf(fid_rp,'Minimum & Maximum limits for model parameter, P_M\n');
fprintf(fid_rp,'%10.4f\n',0.0);
fprintf(fid_rp,'%10.4f\n',0.3);

fprintf(fid_rp,'Minimum & Maximum limits for model parameter, P_Q\n');
fprintf(fid_rp,'%10.4f\n',0.08);
fprintf(fid_rp,'%10.4f\n',0.12);

fprintf(fid_rp,'Minimum & Maximum limits for model parameter, P_SR\n');
fprintf(fid_rp,'%10.4f\n',0.5);
fprintf(fid_rp,'%10.4f\n',1.0);

fclose(fid_rp);


end
