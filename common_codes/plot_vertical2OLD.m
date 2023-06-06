function plot_vertical2OLD(xlims,Important_dates)


    datetick('x','mmm');
    set(gca,'xlim',xlims);
    
    datetick('x','mmm','keeplimits');
    set(gca,'Ygrid','on');
  
    aa=axis;
    for idates=1:length(Important_dates)
        date000=datenum(Important_dates(idates,:));
        plot([date000,date000], [aa(3),aa(4)], 'k-','color',[0.7,0.7,0.7],'LineWidth',1.0);
    end
%     
%     for idates=1:length(Important_dates2)
%         date000=datenum(Important_dates2(idates,:));
%         plot([date000,date000], [aa(3),aa(4)], 'r--');
%     end
%     
%     if(yeslabel==1)
%         for idates=1:length(Important_dates)
%             date000=datenum(Important_dates(idates,:));
%             text(date000,labelY,Important_dates_labels{idates},...
%                 'HorizontalAlignment','left',...
%                 'Rotation',90,...
%                 'VerticalAlignment','bottom',...
%                 'Fontweight','bold',...
%                 'Fontsize',11,...
%                 'Color',[139/255,0,0]);
%         end
%     end
end