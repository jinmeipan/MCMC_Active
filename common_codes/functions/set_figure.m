function titlestr=set_figure(years,width,height)

    fig=figure;set(gcf,'Position',[0,0,width,height]);
    set(fig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
    
    subplot('position',[0.1,0.1,0.8,0.8]);
    
    switch years
        case 2009
            title('IOP1 2009-2010')
            titlestr='IOP1 2009-2010';
        case 2010
            title('IOP2 2010-2011')
            titlestr='IOP2 2010-2011';
        case 2011
            title('IOP3 2011-2012')
            titlestr='IOP3 2011-2012';
        case 2012
            title('IOP4 2012-2013')
            titlestr='IOP4 2012-2013';
    end
end