%plot density,grain size&temperature
function sdkl_plot_snowlayer(years,property)

    switch years
        case 2009
            iop=1;
        case 2010
            iop=2;
        case 2011
            iop=3;
        case 2012
            iop=4;
    end
    
    cd('/Users/jinmeipan/Downloads/MCMC_DataPrep/codes_Read_Snowpit/sodankyla/')
    eval(['load(''sp_sdkl_iop',num2str(iop),'_V0.mat'');']);
    eval(['sps=sp_sdkl_iop',num2str(iop),';']);
    
    %plot snow properties
    for i=1:length(sps)
        sp=sps(i);
        
        dz=sp.dz*100;
        
        %calculate heights
        y(1)=0;
        for i=1:length(dz)
            y(i+1)=y(i)+dz(i);
        end
        
        %
        eval(['val=sp.',property,';']);
        if(strcmp(property,'T')==1)
           val=val-273.15;
           val=max(val,-24.8);
        end
        
        if(sum(isnan(val))==sp.nlayer)
            continue
        end
        
        switch property
            case 'dmax'
                max_D=4;
                colorstr=jet(max_D/0.25);
                %colorstr=flipud(colorstr);
                val_step=0.25;
            case 'density'
                max_dens=500;
                colorstr=jet(max_dens/20);
                colorstr=flipud(colorstr);
                val_step=20;
            case 'T'
                max_T=5;
                min_T=-25;
                NN=((max_T-min_T)/0.2)/3*4;
                jet_temp=jet(NN);
                colorstr=jet_temp([1:125,end-24:end],:);
                val_step=0.2;
        end
        
        %plot layers and fill-in colors
        for i=1:length(dz)
            date0=datenum(sp.year,sp.month,sp.date,sp.time,0,0);
            pos=[date0-1,y(i),2,dz(i)];
            
            if(isnan(val(i))~=1 & val(i)~=0)
                
                ic=round(val(i)/val_step);
                if(strcmp(property,'T')==1)
                    ic=round((val(i)-min_T)/val_step);
                end
                
                rectangle('Position',pos,'FaceColor',colorstr(ic,:),'EdgeColor','k',...
                'LineWidth',1.0);
            else
                rectangle('Position',pos,'FaceColor','w','EdgeColor','k',...
                'LineWidth',1.0);
            end
        end
    end
    
    
    colormap(colorstr);
    h=colorbar('FontSize',12);

    switch property
        case 'dmax'
            set(h, 'ylim', [0 1]);
            set(h, 'YTick',[0:0.5:max_D]/max_D,'YTickLabel',...
            {'0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'});
            set(get(h,'title'),'string','D_{max} (mm)','FontSize',12);
            set(h, 'Position',[0.834,0.285,0.025,0.48])
        case 'density'
            set(h, 'ylim', [0 1]);
            set(h, 'YTick',[0:100:max_dens]/max_dens,'YTickLabel',...
            {'0','100','200','300','400','500'});
            set(get(h,'title'),'string','Density (kg/m^3)','FontSize',12);
            set(h, 'Position',[0.834,0.285,0.025,0.48])   
        case 'T'
            set(h, 'ylim', [0 1]);
            set(h, 'YTick',[0:5:max_T-min_T]/(max_T-min_T),'YTickLabel',...
            {'-25','-20','-15','-10','-5','0','1'});
            set(get(h,'title'),'string','T (^oC)','FontSize',12);
            set(h, 'Position',[0.834,0.285,0.025,0.48])
    end
        
end