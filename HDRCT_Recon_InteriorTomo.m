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


%%
sg = sino_geom('par', 'nb', 1000, 'na', 360, 'dr', 1);
img = LAC_380keV(:,:,ceil(end/2));
ig = image_geom('nx', size(img,1), 'ny', size(img,2), 'fov', size(img,1)*imgdx);

G = Gtomo2_strip(sg, ig);
li = G * img;
figure;imshow(li,[])

img_fbp = em_fbp_QL(sg, ig, li);
figure;imshow(img_fbp,[])

Ind = 1:size(li,1);


gamma = 0;
mu = 1e-05;
[img_TV, Maincost_TV] = IterRecon_TV_FISTA_Simple_2D (Ind, G,li(:), gamma, mu, [ig.nx, ig.ny]);
figure;imshow(img_TV,[])

li_int = li;
Ind = [1:370,630:size(li,1)];
li_int(Ind,:) = 0;
figure;imshow(li_int,[])

img_fbp_int = em_fbp_QL(sg, ig, li_int);
figure;imshow(img_fbp_int,[0,0.15])

gamma = 0;
mu = 1e-05;
[img_int_TV, Maincost_int_TV] = IterRecon_TV_FISTA_Simple_2D (Ind, G,li_int(:), gamma, mu, [ig.nx, ig.ny]);
figure;imshow(img_int_TV,[])


figure;imshow([img_TV img_int_TV],[0,0.2])

figure;imshow([img_fbp, img_fbp_int],[0,0.18])






%% Make measurement - projection

cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250+dx, 'dfs', inf, ...
    'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
    'nt', 1024, 'dt', 500/1024, ...
    'na', 1, 'orbit', 360,'orbit_start',0,...
    'pitch',0, 'source_z0', dz);

ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
    'dx', imgdx, 'dz', imgdz, ...
    'offset_x', dy, 'offset_y', dx, 'offset_z', 0);
figure; cg.plot([ig])

x0 = LAC_380keV;
%         x0(204:207, 51:54, 86) = 10;
%          figure;imshow3D(x0,[])
% figure;imshow(flip(x0(:,:,100),1),[])
f.ptype = 'sf2';
A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
A_cpu.arg.mexfun = @jf_mex;
A_cpus{count} = A_cpu;
y_cpu = A_cpu*x0;
y_cpus{count} = y_cpu;



%% old
img_fbp_noattenuationcorrect = em_fbp_QL(sg, ig, sino);
ind1 = -25; ind2 = 18;
Anni2D = TranslateFigure(Anni3D(:,:,slicenum),ind1,ind2);
% figure;imshow([Anni2D/max(Anni2D(:)),img_fbp_noattenuationcorrect/max(img_fbp_noattenuationcorrect(:))],[])
% C = imfuse(Anni2D/max(Anni2D(:)),img_fbp_noattenuationcorrect/max(img_fbp_noattenuationcorrect(:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% figure; imshow(C)
mumapnew = TranslateFigure(mumap,ind1,ind2);
fluence2new = TranslateFigure(fluence2,ind1,ind2);
maskfluencenew = TranslateFigure(maskfluence,ind1,ind2);

G = Gtomo2_strip(sg, ig);
li = G * mumapnew;
ci = exp(-li);

img_fbp = em_fbp_QL(sg, ig, sino./ci);
img_fbp_corrected = img_fbp./fluence2new;
img_fbp_corrected(maskfluencenew==0)=0;



%
gamma = 0;
mu = 1e-05;

[x_TV, Maincost_TV] = IterRecon_TV_FISTA (x0, A_cpus, y_all(:), gamma, mu, [ig.nx, ig.ny, ig.nz]);



%
gamma = 0;
mu = 1e-05;

[x_TV, Maincost_TV] = IterRecon_TV_FISTA (x0, A_cpus, y_all(:), gamma, mu, [ig.nx, ig.ny, ig.nz]);




%%


% ForBack.applyFP = @(x) for ii = 1:numel(A_cpus) A{ii}*x end;
% ForBack.applyBP = @(x) A'*x;
size1 = A_cpu.size(1);
y_all = zeros(size1*numel(A_cpus),1);
for ii = 1:numel(y_cpus)
    y_all((ii-1)*size1+1:ii*size1) = y_cpus{ii};
end


%
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


