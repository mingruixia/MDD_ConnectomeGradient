function BrainNet_nii2txt(surfacefile,niifile,txtfile)
fid = fopen(surfacefile);
data = textscan(fid,'%f','CommentStyle','#');
surf.vertex_number = data{1}(1);
surf.coord  = reshape(data{1}(2:1+3*surf.vertex_number),[3,surf.vertex_number]);
% ntri = data{1}(3*surf.vertex_number+2);
% tri = reshape(data{1}(3*vertex_number+3:end),[3,ntri])';
fclose(fid);

surf.hdr = spm_vol(niifile);
vol = spm_read_vols(surf.hdr);

surf.coord(4,:)=1;
position = surf.hdr.mat\surf.coord;
position(4,:) = [];
index=round(position);
surf.coord(4,:) = [];

surf.T=zeros(1,surf.vertex_number);
index(:,index(1,:)<=1|index(1,:)>=surf.hdr.dim(1)) = 1;
index(:,index(2,:)<=1|index(2,:)>=surf.hdr.dim(2)) = 1;
index(:,index(3,:)<=1|index(3,:)>=surf.hdr.dim(3)) = 1;

index = sub2ind(surf.hdr.dim,index(1,:),index(2,:),index(3,:));
for i = 1:surf.vertex_number
    if index(i)~=1
        cube = [floor(position(:,i))';ceil(position(:,i))'];
        portion = position(:,i)' - cube(1,:);
        cube(2,portion == 0) = cube(2,portion == 0) + 1;
        tmpT = vol(cube(1,1):cube(2,1),cube(1,2):cube(2,2),cube(1,3):cube(2,3));
        tmpT = (tmpT(:,:,2) - tmpT(:,:,1)) .* portion(3) + tmpT(:,:,1);
        tmpT = (tmpT(:,2) - tmpT(:,1)) .* portion(2) + tmpT(:,1);
        tmpT = (tmpT(2) - tmpT(1)) .* portion(1) + tmpT(1);
        surf.T(i) = tmpT;
    end
end

value = surf.T;
save(txtfile,'value','-ascii');

