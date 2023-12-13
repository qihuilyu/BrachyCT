clear
close all
clc

patname = 'T2';
datatestpath = fullfile('/Users/lyuqihui/Desktop/Data/HDRCT/Datatest/',patname);
load(fullfile(datatestpath,'Datatest.mat'));

imgdx = spatial.PixelSpacings(1);
imgdy = imgdx;
temp = diff(spatial.PatientPositions);
imgdz = temp(1,3);
imgsz = size(CTimage);

offsetpix = 100;
[CenterOfMass, FOV, deltax, deltay, deltaz]  = GetDwellPos(Mask_prostate,4,4,5,imgdx,imgdx,imgdz);

Nfac = 1;
LAC = imresize(LAC_380keV,1/Nfac);
LAC = LAC(:,:,1:Nfac:end);
imgdx = imgdx*Nfac;
imgdy = imgdy*Nfac;
imgdz = imgdz*Nfac;
imgsz = size(LAC);



%% Parameter setting
param.nx = imgsz(1); % number of voxels
param.ny = imgsz(2);
param.nz = imgsz(3);

param.sx = imgdx*param.nx; % mm (real size)
param.sy = imgdy*param.ny; % mm
param.sz = imgdz*param.nz; % mm

%The real detector panel pixel density (number of pixels)
param.nu = 300;		% number of pixels
param.nv = 300;

% Detector setting (real size)
param.su = 500;	% mm (real size)
param.sv = 500;     % mm

% X-ray source and detector setting
param.DSD0 = 250;    %  Distance source to detector
param.DSO0 = 0;	%  X-ray source to object axis distance

% param.DSD = 250;    %  Distance source to detector
% param.DSO = -5;	%  X-ray source to object axis distance

param.dx = imgdx; % single voxel size
param.dy = imgdy;
param.dz = imgdz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_u = 0; param.off_v = 0; % detector rotation shift (real size)

% % % Geometry calculation % % %
param.xs0 = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx;
param.ys0 = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs0 = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz;

param.us0 = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs0 = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;


%% Make measurement - projection
% load img128.mat % Ground-truth image
% M = projection0(param);
% testproj = reshape(M*double(img(:)),[param.nu, param.nv]);
% 
% figure; imshow(testproj,[])

%%

tic

count = 1;
for dz = deltaz
    for dx = deltax
        for dy = deltay

            param.DSD = param.DSD0+dy;    %  Distance source to detector
            param.DSO = param.DSO0+dy;	%  X-ray source to object axis distance

            param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx - dx;
            param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
            param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz - dz;

            param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u - dx;
            param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v - dz;

            M{count} = projection0(param, LAC);
            proj(:,:,count) = reshape(M{count}*double(LAC(:)),[param.nu, param.nv]);

            figure(1);imshow(proj(:,:,count),[])
            % pause(2)
            count = count+1;
        end
    end
end

toc

numproj = count-1;


tic
M_all = [];
for ii = 1:numproj
    M_all = cat(1, M_all, M{ii});
end

proj_all = reshape(M_all*double(LAC(:)),[param.nu, param.nv, numproj]);
figure;imshow3D(proj_all,[])

toc

% figure;gif('Projections.gif')
% for ii = 1:numel(y_cpus)
%     y_all(:,:,ii) = y_cpus{ii};
%     imshow(y_cpus{ii},[12,30]);
%     gif
% end

%% Create projections
x1 = LAC;
y1_all = M_all*x1(:);  
x2 = permute(LAC,[2,1,3]);
y2_all = M_all*x2(:);  
x3 = flip(x1,2);
y3_all = M_all*x3(:);
x4 = flip(x2,2);
y4_all = M_all*x4(:);  

%%

gamma = 0;
mu = 1e-05;


tic
x0 = x1;
x0(x0>0.05) = 0.11;

[x_TV, Maincost_TV] = IterRecon_TV_FISTA_BrachyCT (x0, M_all, y1_all(:), y2_all(:), y3_all(:), y4_all(:), gamma, mu, [param.nx, param.ny, param.nz]);

toc



