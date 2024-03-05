function x_gen_voxelcoord(maskfile,voxelcoordfile)
% maskfile: a binarized nii file of mask
% voxelcoordfile: file to save coordinates for voxels in the mask
mask_hdr = spm_vol(maskfile);
[mask_vol,XYZ] = spm_read_vols(mask_hdr);
coord = XYZ(:,mask_vol>0)';
save(voxelcoordfile,'coord','-ascii');