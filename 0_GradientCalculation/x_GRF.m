function [ClusterSize,dLh,FWHM] = x_GRF(R_volume,DOF,MaskData,Vox,VoxelPThreshold,ClusterPThreshold,tail)
R_Volume = single(R_volume);
[N1, N2, N3, N4]=size(R_Volume);
R_Volume=(R_Volume-repmat(mean(R_Volume,4),[1,1,1, N4]))./repmat(std(R_Volume,0,4),[1,1,1, N4]);%Zero mean and one std
R_Volume(isnan(R_Volume))=0;

SSminus=[0 0 0];
S2=[0 0 0];

N=0;
for x=2:N1
    for y=2:N2
        for z=2:N3
            if MaskData(x, y, z) && MaskData(x-1, y, z) && MaskData(x, y-1, z) && MaskData(x, y, z-1)
                N=N+1;
                for t=1:N4
                    SSminus(1) = SSminus(1) + R_Volume(x, y, z, t) * R_Volume(x-1, y, z, t);
                    SSminus(2) = SSminus(2) + R_Volume(x, y, z, t) * R_Volume(x, y-1, z, t);
                    SSminus(3) = SSminus(3) + R_Volume(x, y, z, t) * R_Volume(x, y, z-1, t);
                    
                    S2(1) = S2(1) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x-1, y, z, t)^2));
                    S2(2) = S2(2) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y-1, z, t)^2));
                    S2(3) = S2(3) + 0.5 * ((R_Volume(x, y, z, t)^2) + (R_Volume(x, y, z-1, t)^2));
                end
            end
        end
    end
    %     fprintf('.');
end

if SSminus(1)>0.99999999*S2(1)
    SSminus(1)=0.99999999*S2(1);
    warning('possibly biased smootheness in X');
end
if SSminus(2)>0.99999999*S2(2)
    SSminus(2)=0.99999999*S2(2);
    warning('possibly biased smootheness in Y');
end
if SSminus(3)>0.99999999*S2(3)
    SSminus(3)=0.99999999*S2(3);
    warning('possibly biased smootheness in Z');
end

sigmasq(1) = -1 / (4 * log(abs(SSminus(1)/S2(1))));
sigmasq(2) = -1 / (4 * log(abs(SSminus(2)/S2(2))));
sigmasq(3) = -1 / (4 * log(abs(SSminus(3)/S2(3))));

dLh=((sigmasq(1)*sigmasq(2)*sigmasq(3))^-0.5)*(8^-0.5);

if N4 > 1
%     fprintf('DLH %f voxels^-3 before correcting for temporal DOF\n',dLh);
    
    lut(6)   = 1.5423138; lut(7)   = 1.3757105; lut(8)   = 1.2842680;
    lut(9)   = 1.2272151; lut(10)  = 1.1885232; lut(11)  = 1.1606988;
    lut(12)  = 1.1398000; lut(13)  = 1.1235677; lut(14)  = 1.1106196;
    lut(15)  = 1.1000651; lut(16)  = 1.0913060; lut(17)  = 1.0839261;
    lut(18)  = 1.0776276; lut(19)  = 1.0721920; lut(20)  = 1.0674553;
    lut(21)  = 1.0632924; lut(26)  = 1.0483053; lut(31)  = 1.0390117;
    lut(41)  = 1.0281339; lut(51)  = 1.0219834; lut(61)  = 1.0180339;
    lut(71)  = 1.0152850; lut(81)  = 1.0132621; lut(91)  = 1.0117115;
    lut(101) = 1.0104851; lut(151) = 1.0068808; lut(201) = 1.0051200;
    lut(301) = 1.0033865; lut(501) = 1.0020191;
    
    y = lut(lut~=0);
    x = find(lut~=0);
    xi=[1:501];
    lut_interpolated=interp1(x,y,xi,'linear');
    
    if (DOF < 6)
        dLh=dLh * 1.1;
    elseif (DOF>500)
        dLh=dLh * (1.0321/DOF +1)^0.5;
    else
        retval=(lut_interpolated(floor(DOF)+1)-lut_interpolated(floor(DOF)))*(floor(DOF)+1-floor(DOF)) + ...
            lut_interpolated(floor(DOF)+1);
        dLh=dLh * retval^0.5;
    end
    
end

FWHM(1) =  sqrt(8 * log(2) * sigmasq(1));
FWHM(2) =  sqrt(8 * log(2) * sigmasq(2));
FWHM(3) =  sqrt(8 * log(2) * sigmasq(3));

resels = FWHM(1)*FWHM(2)*FWHM(3);
% fprintf('\nFWHMx = %f voxels\nFWHMy = %f voxels\nFWHMz = %f voxels\n',FWHM(1),FWHM(2),FWHM(3));
FWHM=FWHM.*Vox;
% fprintf('FWHMx = %f mm\nFWHMy = %f mm\nFWHMz = %f mm\n',FWHM(1),FWHM(2),FWHM(3));
nVoxels=length(find(MaskData));
% fprintf('DLH = %f\nVOLUME = %d\nRESELS = %f\n',dLh,nVoxels,resels);


if tail == 2
    ClusterPThreshold = ClusterPThreshold/2;
    VoxelPThreshold = VoxelPThreshold/2;
end


zThrd=norminv(1 - VoxelPThreshold);

D=3;
Em = nVoxels * (2*pi)^(-(D+1)/2) * dLh * (zThrd*zThrd-1)^((D-1)/2) * exp(-zThrd*zThrd/2);
EN = nVoxels * (1-normcdf(zThrd)); %In Friston et al., 1994, EN = S*Phi(-u). (Expectation of N voxels)  % K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.
Beta = ((gamma(D/2+1)*Em)/(EN)) ^ (2/D); % K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.

% Get the minimum cluster size
pTemp=1;
ClusterSize=0;
while pTemp >= ClusterPThreshold
    ClusterSize=ClusterSize+1;
    pTemp = 1 - exp(-Em * exp(-Beta * ClusterSize^(2/D))); %K. Friston, K. Worsley, R. Frackowiak, J. Mazziotta, and A. Evans. Assessing the significance of focal activations using their spatial extent. Human Brain Mapping, 1:214?220, 1994.
end

% fprintf('%f voxels\n',ClusterSize);

