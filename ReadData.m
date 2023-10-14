close all
clc
clear

patname = 'T2';
datapath = fullfile('/Users/qihuilyu/Desktop/Data/HDRCT/RawData/', patname);
ctpath = fullfile(datapath,'CT');
RTdosepath = fullfile(datapath,'RTDOSE');
RTplanpath = fullfile(datapath,'RTPLAN');
RTstpath = fullfile(datapath,'RTst');
datatestpath = fullfile('/Users/qihuilyu/Desktop/Data/HDRCT/Datatest/',patname);
mkdir(datatestpath)

[CTimage,spatial,dim] = dicomreadVolume(ctpath);
CTfiles = dir(fullfile(ctpath,'**/*.dcm'));
info = dicominfo(fullfile(ctpath,CTfiles(1).name));
CTimage = (flip(squeeze(double(CTimage)),1)+info.RescaleIntercept)*info.RescaleSlope;
figure;imshow3D(CTimage,[-500,500])

RTstfile = dir(fullfile(RTstpath,'**/*.dcm'));
file = fullfile(RTstpath, RTstfile.name);
infoRTst = dicominfo(file,'UseVRHeuristic',false);
rtContours = dicomContours(infoRTst);
plotContour(rtContours)

% referenceInfo = imref3d(size(V),xlim,ylim,zlim);
for ii = 1:size(rtContours.ROIs,1)
    rtMasks{ii,1} = flip(createMask(rtContours, ii, spatial),1);
    StructureVoxelNum{ii,1} = nnz(rtMasks{ii,1});
end

StructureName = rtContours.ROIs.Name;
StructureInfo = struct('Name',StructureName,'Mask',rtMasks,...
    'VoxelNum',StructureVoxelNum);


Mask_prostate = StructureInfo(13).Mask;
figure;imshow3D(Mask_prostate,[])

imgdx = spatial.PixelSpacings(1);
temp = diff(spatial.PatientPositions);
imgdz = temp(1,3);
[CenterOfMass, FOV] = GetPTV_COM_FOV(Mask_prostate,imgdx,imgdx,imgdz);

shiftxyz = CenterOfMass-ceil(spatial.ImageSize/2);
gap = [40,40];

NewCT = ones(spatial.ImageSize + [0,0,2*abs(shiftxyz(3))])*(-1024);
NewCT(gap(1)-shiftxyz(1):end-shiftxyz(1)-gap(1),gap(2)-shiftxyz(2):end-shiftxyz(2)-gap(2),abs(shiftxyz(3))+1-shiftxyz(3):end-abs(shiftxyz(3))-shiftxyz(3)) = CTimage(gap(1):end-gap(1),gap(2):end-gap(2),:);

for ii = 1:size(rtContours.ROIs,1)
    OldMask = StructureInfo(ii).Mask;
    NewMask = zeros(spatial.ImageSize + [0,0,2*abs(shiftxyz(3))]);
    NewMask(gap(1)-shiftxyz(1):end-shiftxyz(1)-gap(1),gap(2)-shiftxyz(2):end-shiftxyz(2)-gap(2),abs(shiftxyz(3))+1-shiftxyz(3):end-abs(shiftxyz(3))-shiftxyz(3)) = OldMask(gap(1):end-gap(1),gap(2):end-gap(2),:);

    StructureInfo(ii).Mask = NewMask;
    StructureInfo(ii).VoxelNum = nnz(NewMask);

end

Mask_prostate = StructureInfo(13).Mask;

% Mask_prostateNew = zeros(spatial.ImageSize + [0,0,2*abs(shiftxyz(3))]);
% Mask_prostateNew(gap(1)-shiftxyz(1):end-shiftxyz(1)-gap(1),gap(2)-shiftxyz(2):end-shiftxyz(2)-gap(2),abs(shiftxyz(3))+1-shiftxyz(3):end-abs(shiftxyz(3))-shiftxyz(3)) = Mask_prostate(gap(1):end-gap(1),gap(2):end-gap(2),:);
% [CenterOfMass, FOV] = GetPTV_COM_FOV(Mask_prostateNew,imgdx,imgdx,imgdz);
%


% figure;imshow3D([rtMasks{21}, double(CTimage)/double(max(CTimage(:)))+1],[]);

CTimage = NewCT;

save(fullfile(datatestpath,'Datatest.mat'),'CTimage','StructureInfo','spatial','Mask_prostate');


