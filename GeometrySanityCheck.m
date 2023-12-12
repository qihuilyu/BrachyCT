addpath('bin');

ParamSetting;

%% Make measurement - projection
load img128.mat % Ground-truth image

tic
proj = CTprojection(img,param);
toc
proj3D = reshape(proj, [256,200,360]); 
figure; imshow3D(proj3D,[0,6])

save proj.mat proj 

%% Sanity check
cg = ct_geom('fan', 'dsd', 1500, 'dod', 400, 'dfs', inf, ...
    'ns', 256, 'ds', 4, 'offset_s', 0, ...
    'nt', 200, 'dt', 4, ...
    'na', 360, 'orbit', 360,'orbit_start',0,...
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

figure; imshow3D([proj3D flip(proj2,1)],[0,6])
figure; imshow3D([proj3D-flip(proj2,1)],[0,6])


%% Check A': <Ax,y> = <A'y,x>

x = rand(size(img));
x([1:25,104:128],:,:) = 0;
x(:,[1:25,104:128],:) = 0;
x(:,:,[1:25,104:128]) = 0;


Ax = CTprojection(x,param);
y = rand(size(proj));
Aty = CTbackprojection(y, param);

% Norimg = CTbackprojection(CTprojection(ones(param.nx,param.ny,param.nz,'single'),param), param);
(sum(Ax(:).*y(:)) - sum(Aty(:).*x(:)))/sum(Ax(:).*y(:))
(sum(Ax(:).*y(:)) - sum(Aty(:).*x(:)))


x2 = flip(x,1);
Ax2 = A_cpu*x2;
y2 = flip(y,1);
Aty2 = A_cpu'*y2;
sum(Ax2(:).*y2(:)) - sum(Aty2(:).*x2(:)) 




test1 = proj3D-flip(proj2,1); norm(test1(:))/norm(proj3D(:))

test1 = Ax-flip(Ax2,1); norm(test1(:))/norm(Ax(:))
figure;imshow3D([Ax, flip(Ax2,1),Ax-flip(Ax2,1)],[])
figure;imshow3D(Ax-flip(Ax2,1),[])

test1 = Aty2-flip(Aty,1); norm(test1(:))/norm(Aty2(:))
figure;imshow3D([Aty2, flip(Aty,1)*5,Aty2-flip(Aty,1)*5],[])
figure;imshow3D(Aty2-flip(Aty,1)*5,[])

