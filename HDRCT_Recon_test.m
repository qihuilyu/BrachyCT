clear
close all




dz = 0;

for dx = [-100,0,100]


cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, ...
    'ns', 1024, 'ds', 400/1024, 'offset_s', 0, ...
    'nt', 1024, 'dt', 400/1024, ...
    'na', 1, 'orbit', 360,'orbit_start',0,...
    'pitch',0, 'source_z0', dz);

ig = image_geom('nx', 512, 'ny', 512, 'nz', 200, ...
    'dx', 1, 'dz', 1, ...
    'offset_x', 0, 'offset_y', 0, 'offset_z', 0);

figure; cg.plot([ig])	

x0 = zeros(ig.nx, ig.ny, ig.nz);
x0(180:210,250:262,:) = 1;
figure;imshow(flip(x0(:,:,100),1),[])
f.ptype = 'sf2';
A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
A_cpu.arg.mexfun = @jf_mex;
y_cpu = A_cpu*x0;
figure;imshow(y_cpu,[])



end





% 
% 
% 
% f.ptype = 'sf2';
% A_gpu = Gcone_gpu(cg, ig, 'type', f.ptype, 'ngpu', 4);
% A_gpu.arg.mexfun = @cbct_mex;
% 




% cpu etic
% y_gpu = A_gpu*x0;
% cpu etoc 'new projection gpu time: '
% 
% nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
% printm('nrmse proj %g %%', nrmse(y_cpu, y_gpu))

%proj = cg.zeros;
%for ia = 1:cg.na
%    proj(mod(ia,cg.ns)+1, mod(ia,cg.nt)+1, ia) = ia/10;
%end
% x_cpu = A_cpu' * y_cpu;
% cpu etic
% x_gpu = A_gpu' * y_gpu;
% cpu etoc 'back projection gpu time: '
% 
% %temp1 = x_cpu(:);
% 
% %for i = 87000000:87001000
% %    if temp1(i,:) ~= 0.0
% %	printm('%d,%g\n', i, temp1(i,:));
% %    end
% %end
% 
% printm('nrmse back %g %%', nrmse(x_cpu, x_gpu));
