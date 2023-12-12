addpath('bin');

ParamSetting;

param.DSD = 250;    %  Distance source to detector 
param.DSO = -5;	%  X-ray source to object axis distance
param.dang = 45; % angular step size (deg)
param.deg = 0:param.dang:359; % you can change
param.deg = param.deg*param.dir;
param.nProj = length(param.deg);

%% Make measurement - projection
load img128.mat % Ground-truth image

tic
proj = CTprojection(img,param);
toc
proj3D = reshape(proj, [256,200,8]); 
figure; imshow3D(proj3D,[0,6])

save proj.mat proj 

%% Sanity check
% It seems that the IRT toolbox is not useful in BrachyCT geometry (when the source is inside the object)
% Under standard geometry the matlab projection toolbox matches with the
% irt toolbox
cg = ct_geom('fan', 'dsd', param.DSD, 'dod', param.DSD - param.DSO, 'dfs', inf, ...
    'ns', 256, 'ds', 4, 'offset_s', 0, ...
    'nt', 200, 'dt', 4, ...
    'na', 8, 'orbit', 360,'orbit_start',0,...
    'pitch',0, 'source_z0', 0);

ig = image_geom('nx', 128, 'ny', 128, 'nz', 128, ...
    'dx', 460/128, 'dz', 460/128, ...
    'offset_x', 0, 'offset_y', 0, 'offset_z', 0);
figure; cg.plot([ig])

f.ptype = 'sf2';

A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
A_cpu.arg.mexfun = @jf_mex;

tic
proj2 = A_cpu*flip(img,1);
toc

figure; imshow3D([proj3D flip(proj2,1)],[0,10])
figure; imshow3D([proj3D-flip(proj2,1)],[-6,6])
