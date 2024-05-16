clear
close all
clc
profile on

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
NfacZ = 1;
LAC = imresize(LAC_380keV,1/Nfac);
LAC = LAC(:,:,1:NfacZ:end);
imgdx = imgdx*Nfac;
imgdy = imgdy*Nfac;
imgdz = imgdz*NfacZ;
imgsz = size(LAC);



%% Parameter setting
param.nx = imgsz(1); % number of voxels
param.ny = imgsz(2);
param.nz = imgsz(3);

param.sx = imgdx*param.nx; % mm (real size)
param.sy = imgdy*param.ny; % mm
param.sz = imgdz*param.nz; % mm

%The real detector panel pixel density (number of pixels)
param.nu = 1000;		% number of pixels
param.nv = 1000;

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
dz = 0;
dx = 0;
dy = 0;

param.DSD = param.DSD0+dy;    %  Distance source to detector
param.DSO = param.DSO0+dy;	%  X-ray source to object axis distance

param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx - dx;
param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz - dz;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u - dx;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v - dz;

% M{count} = projection0(param, LAC);
% proj(:,:,count) = reshape(M{count}*double(LAC(:)),[param.nu, param.nv]);

proj(:,:) = projection0_Nomatrix(param, LAC);

figure(1);imshow(proj(:,:),[])


%%
LAC_new = LAC;

BBradius = 2;

BBcenter = [400,300,ceil(size(LAC,3)/2)];
[ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
LAC_new(ind) = 3.86;

BBcenter = [403,200,ceil(size(LAC,3)/2)];
[ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
LAC_new(ind) = 3.86;


BBradius = 2.8;

BBcenter = [425,350,ceil(size(LAC,3)/2)-20];
[ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
LAC_new(ind) = 3.86;

BBcenter = [425,250,ceil(size(LAC,3)/2)-20];
[ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
LAC_new(ind) = 3.86;

% BBcenter = [400,300,ceil(size(LAC,3)/2)-20];
% [ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
% ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
% LAC_new(ind) = 3.86;

LAC_new = permute(LAC_new,[2,1,3]);
figure(2);imshow3D(LAC_new,[0,0.15])
figure(10);imshow3D(permute(LAC_new,[2,3,1]),[0,0.2])

proj1(:,:) = projection0_Nomatrix(param, LAC_new);

% figure(3);imshow(proj1(:,:),[])

[centers,radii] = imfindcircles(proj1,[2 10]);

figure; imshow(proj1,[14,23])
h = viscircles(centers,radii);


%%


LAC_new = LAC;


BBradius = 8;
BBcenter = [300,260,ceil(size(LAC,3)/2)+7];
[ind1,ind2,ind3] = ndgrid(1:size(LAC,1),1:size(LAC,2),1:size(LAC,3));
ind = find((ind1-BBcenter(1)).^2*imgdx^2+(ind2-BBcenter(2)).^2*imgdy^2+(ind3-BBcenter(3)).^2*imgdz^2<BBradius^2);
LAC_new(ind) = 0.02;
figure(2);imshow3D([LAC, LAC_new],[0,0.15])


% LAC_new = permute(LAC_new,[2,1,3]);
% figure(10);imshow3D(permute(LAC_new,[2,3,1]),[0,0.2])

% proj1(:,:) = projection0_Nomatrix(param, LAC_new);
% proj0(:,:) = projection0_Nomatrix(param, LAC);
% 

proj1(:,:) = projection0_Nomatrix(param, permute(LAC_new,[2,1,3]));
proj0(:,:) = projection0_Nomatrix(param, permute(LAC,[2,1,3]));


figure(3);imshow([proj0,proj1],[])


%%
ind = 750; LineWidth = 3;
FigureColor = 'white';
set(gcf,'color',FigureColor);
LegendSize = 13;
LabelSize = 15;
LegendLocation = 'northeast';
AxesSize = 15;
TitleSize = 15;


figure;

plot(proj0(:,ind),"--",'LineWidth',LineWidth); hold on;
plot(proj1(:,ind),":",'LineWidth',LineWidth);

    h=legend({'Baseline', 'New air bubble'},'FontSize',LegendSize);
    set(h,'Location',LegendLocation);
    xlabel('Detector pixel','FontSize',LabelSize);
    ylabel('Line Integral','FontSize',LabelSize);
    hline = findobj(gcf, 'type', 'line');
    set(hline,'LineWidth',LineWidth)
    set(gca,'fontsize',AxesSize);
    title('Detector signals with and without simulated air bubbles','FontSize',TitleSize);





%%



toc


% 
% tic
% M_all = [];
% for ii = 1:numproj
%     M_all = cat(1, M_all, M{ii});
% end
% 
% proj_all = reshape(M_all*double(LAC(:)),[param.nu, param.nv, numproj]);
% figure;imshow3D(proj_all,[])
% 
% toc
% 
% % figure;gif('Projections.gif')
% % for ii = 1:numel(y_cpus)
% %     y_all(:,:,ii) = y_cpus{ii};
% %     imshow(y_cpus{ii},[12,30]);
% %     gif
% % end
% 
% %% Create projections
% x1 = LAC;
% y1_all = M_all*x1(:);
% x2 = permute(LAC,[2,1,3]);
% y2_all = M_all*x2(:);
% x3 = flip(x1,2);
% y3_all = M_all*x3(:);
% x4 = flip(x2,2);
% y4_all = M_all*x4(:);
% 
% %%
% 
% gamma = 0;
% mu = 1e-05;
% 
% 
% tic
% x0 = x1;
% x0(x0>0.05) = 0.11;
% 
% [x_TV, Maincost_TV] = IterRecon_TV_FISTA_BrachyCT (x0, M_all, y1_all(:), y2_all(:), y3_all(:), y4_all(:), gamma, mu, [param.nx, param.ny, param.nz]);
% 
% toc
% 
% 
% profile report
