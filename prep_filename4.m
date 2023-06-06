%filename_prepare
%
% input:
% -file: the filename for site,month,model-specific choice

function MCMCRun4=prep_filename4(MCMCRun4)

folder0=MCMCRun4.folder;

fid_fn=fopen([folder0,'FILENAME.txt'],'w');

fprintf(fid_fn,['RunParams.txt\n']);
fprintf(fid_fn,['tb_obs.txt\n']);
fprintf(fid_fn,['hyperpar.txt\n']);
fprintf(fid_fn,['acceptance.txt\n']);
fprintf(fid_fn,['post_tb.out\n']);
fprintf(fid_fn,['post_theta.out']);
fclose(fid_fn);

clear fid_fn

end



