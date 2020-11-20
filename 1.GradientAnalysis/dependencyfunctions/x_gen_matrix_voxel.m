function M = x_gen_matrix_voxel(maskfile,funcfile)
mask_hdr = spm_vol(maskfile);
mask_vol = spm_read_vols(mask_hdr);
mask_ind = reshape(mask_vol>0,1,[]);

Nii = nifti(funcfile);
volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
maskCourse = volCourse(:,mask_ind);

M = corr(maskCourse);

