addpath('bin');

ParamSetting;

%% Recon case 1 - Analytic reconstruction: filtered backprojection
% filter='ram-lak','shepp-logan','cosine', 'hamming', 'hann' : (ramp + additional filter)
param.filter='hann'; 
load proj.mat

proj_filtered = filtering(proj,param);
Reconimg = CTbackprojection(proj_filtered, param);

for i=1:param.nz
    figure(2); imagesc(max(Reconimg(:,:,i),0)); axis off; axis equal; colormap gray; colorbar;
    title(num2str(i));
    pause(0.01);
end



figure;imshow3D([t,tsqrt,t_dist/3.7,img],[])
figure;imshow3D([t-tsqrt,t-img, tsqrt-img],[])
figure;imshow3D([t-tsqrt],[])
figure;imshow3D([t-tsqrt],[-0.001,0.001])
param.zs
figure;imshow3D([t-tsqrt,t-img, tsqrt-img],[])
figure;imshow3D([t-tsqrt,t-img, tsqrt-img],[-0.001,0.001])
figure;imshow3D([t-tsqrt,t-img, tsqrt-img],[-0.001,0.001]); colorbar
figure;imshow3D(abs([t-tsqrt,t-img, tsqrt-img, t_dist-img]),[0,0.001]); colorbar


