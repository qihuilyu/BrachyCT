clear
close all
clc

patname = 'T2';
datatestpath = fullfile('/Users/lyuqihui/Desktop/Data/HDRCT/Datatest/',patname);
load(fullfile(datatestpath,'Datatest.mat'));

imgdx = spatial.PixelSpacings(1);
temp = diff(spatial.PatientPositions);
imgdz = temp(1,3);
imgsz = size(CTimage);

offsetpix = 100;

[CenterOfMass, FOV, deltax, deltay, deltaz]  = GetDwellPos(Mask_prostate,4,4,5,imgdx,imgdx,imgdz);


%% Parameter setting
param.nx = 128; % number of voxels
param.ny = 128;
param.nz = 128;

param.sx = 460; % mm (real size)
param.sy = 460; % mm
param.sz = 460; % mm

%The real detector panel pixel density (number of pixels)
param.nu = 256;		% number of pixels
param.nv = 200;

% Detector setting (real size)
param.su = 500;	% mm (real size)
param.sv = 500;     % mm

% X-ray source and detector setting
param.DSD0 = 1500;    %  Distance source to detector
param.DSO0 = 1100;	%  X-ray source to object axis distance

% param.DSD = 250;    %  Distance source to detector
% param.DSO = -5;	%  X-ray source to object axis distance

param.dx = param.sx/param.nx; % single voxel size
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
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
load img128.mat % Ground-truth image
% M = projection0(param);
% testproj = reshape(M*double(img(:)),[param.nu, param.nv]);
% 
% figure; imshow(testproj,[])


%% Sanity check
% tic
% param.dir = -1;   % gantry rotating direction (clock wise/ counter clockwise)
% param.dang = 360; % angular step size (deg)
% param.deg = 0:param.dang:359; % you can change
% param.deg = param.deg*param.dir;
% param.nProj = length(param.deg);
% param.filter='ram-lak'; % high pass filter
% param.interptype = 'linear'; % 'linear', 'nearest'
% param.gpu = 1;
%
% proj2 = CTprojection(img,param);
% figure; imshow([testproj proj2],[])
% toc
% proj3D = reshape(proj, [256,200,8]);
% figure; imshow3D(proj3D,[0,6])
% save proj.mat proj


%%
count = 1;
for dz = [-200,0,200]
    for dx = [-200,0,200]
        for dy = [-200,0,200]

            param.DSD = param.DSD0+dy;    %  Distance source to detector
            param.DSO = param.DSO0+dy;	%  X-ray source to object axis distance

            param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx - dx;
            param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
            param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz - dz;

            param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u - dx;
            param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v - dz;

            M{count} = projection0(param);
            proj(:,:,count) = reshape(M{count}*double(img(:)),[param.nu, param.nv]);

            figure(1);imshow(proj(:,:,count),[])
            pause(2)
            count = count+1;
        end
    end
end

numproj = count-1;

M_all = [];
for ii = 1:numproj
    M_all = cat(1, M_all, M{ii});
end

proj_all = reshape(M_all*double(img(:)),[param.nu, param.nv, numproj]);
figure;imshow3D(proj_all,[])


% figure;gif('Projections.gif')
% for ii = 1:numel(y_cpus)
%     y_all(:,:,ii) = y_cpus{ii};
%     imshow(y_cpus{ii},[12,30]);
%     gif
% end


%%


% ForBack.applyFP = @(x) for ii = 1:numel(A_cpus) A{ii}*x end;
% ForBack.applyBP = @(x) A'*x;
size1 = A_cpu.size(1);
y_all = zeros(size1*numel(A_cpus),1);
for ii = 1:numel(y_cpus)
    y_all((ii-1)*size1+1:ii*size1) = y_cpus{ii};
end

gamma = 0;
mu = 1e-05;

[x_TV, Maincost_TV] = IterRecon_TV_FISTA (x0, A_cpus, y_all(:), gamma, mu, [ig.nx, ig.ny, ig.nz]);

figure;
test = [LAC_380keV x_TV];
for ii = 1:size(test,3)
    imshow(test(:,:,ii),[]);
    exportgraphics(gcf,'testx0.gif','Append',true);
end

figure;
test = [reshape(Ax0,[1024,1024,32]) reshape(y_all,[1024,1024,32])];
for ii = 1:size(test,3)
    imshow(test(:,:,ii),[]);
    exportgraphics(gcf,'test.gif','Append',true);
end



% x0 = zeros(ig.nx, ig.ny, ig.nz);
% x0(180:210,250:262,:) = 1;
% figure;imshow(flip(x0(:,:,100),1),[])
% f.ptype = 'sf2';
% A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
% A_cpu.arg.mexfun = @jf_mex;
% y_cpu = A_cpu*x0;
% figure;imshow(y_cpu,[])


