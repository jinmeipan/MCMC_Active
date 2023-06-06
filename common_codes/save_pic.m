
% function save_pic(filename)
% save it on mac desktop in jpg format
% Input:

% filename=string, eg, 'fig1'
% '.jpg' will be added automatically

function save_pic(filename,folder)


set(gcf, 'PaperPositionMode', 'auto')
set(gcf,'color','w');

if(nargin<2)
% folder='/Users/jinmeipan/Desktop/Simu_Sodankyla/IOPs_5/';
% folder='/Users/jinmeipan/Downloads/MCMCRunData_V3/R1/2lyr-syns_0.1/';
% folder='/Users/jinmeipan/Desktop/Simu_Sodankyla/IOPs_3/new figures/';
% folder='/Users/jinmeipan/Downloads/MCMCRunData_V3/R1/thesis figures/';
% folder='/Users/jinmeipan/Downloads/MCMCRunData_V3/R1/!sysn_1lyr_res/';
% folder='/Users/jinmeipan/Downloads/MCMCRunData_V3_PMtest/Longmont_pic/';
folder='/Users/jinmeipan/Desktop/untitled folder/';
end
mkdir(folder);

name=[folder,filename,'.eps']
name2=[folder,filename,'.jpg'];

% % print('-depsc','-painters',name)
% % print('-djpeg','-r300',name2)


savefig([folder,'\',filename,'.fig']);
print('-djpeg','-r600',[folder,'\',filename,'.jpg']);
print('-dtiff','-r800',[folder,'\',filename,'.tif']);
print('-depsc',[folder,'\',filename,'.eps']);


end
