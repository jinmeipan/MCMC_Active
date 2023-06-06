%inputs:
%tb_type0: 1=real Tb, 2=simulated Tb
%note, the variable TB_MCMC actually allows to contain both active and
%passive microwave observations


function MCMCRun4=prep_tbobs4(MCMCRun4,tb_type0,PitNo)


folder0=MCMCRun4.folder;


switch tb_type0
    case 'real'
        tb_type=1;
        disp('Used real Tb')
    otherwise
        tb_type=2;
        disp('Used synthetic Tb')
end

freq=MCMCRun4.MCMC_freq; nfreq=length(freq);
theta=MCMCRun4.MCMC_theta; ntheta=length(theta);
np=MCMCRun4.passive_np + MCMCRun4.active_np + MCMCRun4.other_np;

%read or caculate TB_MCMC
%note, don't input data which is not given for freq&theta for each specific
%sensor
sp=MCMCRun4.sp;

    
%------

TB_MCMC0=zeros(np,ntheta,nfreq)+999; %default value, 999

for iff=1:nfreq
    for itt=1:ntheta
        for ipp=1:np
            %add data, for active
            if(ipp>0 & ipp<=MCMCRun4.passive_np)
                %it is in passive
                idx1=find(freq(iff)==MCMCRun4.passive_freq);
                idx2=find(theta(itt)==MCMCRun4.passive_theta);
                if(length(idx1)==0 | length(idx2)==0)
                    continue
                end
                
                idx4=find(freq(iff)==sp.freq);
                idx3=find(theta(itt)==sp.theta);
                
                if(tb_type==1)
                    temptb=[sp.tbv(idx3,idx4),sp.tbh(idx3,idx4)];
                else
                    temptb=[sp.tbv_simu(idx3,idx4),sp.tbh_simu(idx3,idx4)];
                    temptb=temptb+random('Normal',0,MCMCRun4.stdTb,size(temptb,1),size(temptb,2));
                end
                    
                switch MCMCRun4.passive_pol{ipp}
                    case 'v'
                        TB_MCMC0(ipp,itt,iff)=temptb(1);
                    case 'h'
                        TB_MCMC0(ipp,itt,iff)=temptb(2);
                end
            end
            
            %%
            if(ipp>MCMCRun4.passive_np & ipp<=MCMCRun4.passive_np+MCMCRun4.active_np)
                %it is in active
                idx1=find(freq(iff)==MCMCRun4.active_freq);
                idx2=find(theta(itt)==MCMCRun4.active_theta);
                if(length(idx1)==0 | length(idx2)==0)
                    continue
                end
                
                idx4=find(freq(iff)==sp.freq_ac);
                idx3=find(theta(itt)==sp.theta_ac);
                
                if(tb_type==1)
                    temptb=[sp.sigma_vv(idx3,idx4),sp.sigma_hh(idx3,idx4),...
                            sp.sigma_vh(idx3,idx4),sp.sigma_hv(idx3,idx4)];
                else
                    temptb=[sp.sigma_vv_simu(idx3,idx4),sp.sigma_hh_simu(idx3,idx4),...
                            sp.sigma_vh_simu(idx3,idx4),sp.sigma_hv_simu(idx3,idx4)];
                    temptb=temptb+random('Normal',0,MCMCRun4.stdSigma,size(temptb,1),size(temptb,2));
                end
                
                switch MCMCRun4.active_pol{ipp+MCMCRun4.passive_np}
                    case 'vv'
                        TB_MCMC0(ipp,itt,iff)=temptb(1);
                    case 'hh'
                        TB_MCMC0(ipp,itt,iff)=temptb(2);
                    case 'vh'
                        TB_MCMC0(ipp,itt,iff)=temptb(3);
                    case 'hv'
                        TB_MCMC0(ipp,itt,iff)=temptb(4);
                end
            end
            
            
            %%
            if(ipp>MCMCRun4.passive_np+MCMCRun4.active_np & iff==1)
                TB_MCMC0(ipp,itt,iff)=MCMCRun4.other_measurements;
            end
        end
    end
end


TB_MCMC=TB_MCMC0(:,1,1)';
for iff=1:nfreq
    for itt=1:ntheta
        if(iff==1 & itt==1)
            continue
        end
        
        TB_MCMC=[TB_MCMC;TB_MCMC0(:,itt,iff)'];
    end
end

MCMC4.ObsMCMCIn=TB_MCMC;


%% write
fid_tb=fopen([folder0,'tb_obs.txt'],'w');
fprintf(fid_tb,[' Pit  #     ',1,'  \n']);

for i=1:size(TB_MCMC,1)
    
    strm=repmat('%12.4f',1,np);
    strm=[strm,'\n'];
    
    fprintf(fid_tb,strm, TB_MCMC(i,:));
end
fclose(fid_tb);


end
