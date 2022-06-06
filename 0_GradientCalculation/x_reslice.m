function x_reslice(refimg,sourceimg,interp)
spmjob = load('Coregister_Reslice.mat');
spmjob.matlabbatch{1,1}.spm.spatial.coreg.write.ref = {[refimg,',1']};
hdr = spm_vol(sourceimg);
spmjob.matlabbatch{1,1}.spm.spatial.coreg.write.source = cell(length(hdr),1);
for i = 1:length(hdr)
    spmjob.matlabbatch{1,1}.spm.spatial.coreg.write.source{i} = [sourceimg,',',num2str(i)];
end
spmjob.matlabbatch{1,1}.spm.spatial.coreg.write.roptions.interp = interp;
spm_jobman('initcfg');
spm_jobman('run',spmjob.matlabbatch);
