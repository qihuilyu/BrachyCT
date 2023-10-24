clear
close all
clc

patname = 'T2';
datatestpath = fullfile('/Users/lyuqihui/Desktop/Data/HDRCT/Datatest/',patname);
load(fullfile(datatestpath,'Datatest.mat'));

% Materials = [Air, Dry
% Adipose Tissue
% Breast Tissue
% water
% Muscle, Skeletal
% Lung Tissue
% Tissue, Soft
% bone cortical]
%

MU_40 = [0.00030
0.22762
0.26830
0.28193
0.28340
0.28493
1.27776];

LAC_40keV = (CTimage/1000+1)*MU_40(3);
LAC_40keV(LAC_40keV<0) = 0;

MU_380 = [0.00012
0.10327
0.10860
0.11306
0.11317
0.11412
0.19493];

LAC_380keV = interp1(MU_40,MU_380,LAC_40keV,'linear','extrap');

save(fullfile(datatestpath,'Datatest.mat'),'LAC_40keV','LAC_380keV',"-append")

figure;imshow3D([LAC_40keV LAC_380keV],[0,0.7])
figure;imshow3D([LAC_380keV],[0.05,0.15])
figure;imshow3D([LAC_40keV],[0,0.7])

%% Sanity check

testslice = 64;
LAC_40keV2D = squeeze(LAC_40keV(:,:,testslice));
[minV_40keV2D, minInd_40keV2D] = min(abs(LAC_40keV2D - reshape(MU_40,[1,1,numel(MU_40)])),[],3);

LAC_380keV2D = squeeze(LAC_380keV(:,:,testslice));
[minV_380keV2D, minInd_380keV2D] = min(abs(LAC_380keV2D - reshape(MU_380,[1,1,numel(MU_380)])),[],3);


figure;imshow([minInd_40keV2D],[])
figure;imshow([minInd_380keV2D],[])
figure;imshow([minInd_40keV2D, minInd_380keV2D],[])
figure;imshow([minInd_40keV2D-minInd_380keV2D],[])




