function proj2d = projection0(data3d,param)

nu = param.nu;
nv = param.nv;

proj2d = (zeros(param.nu,param.nv,'single'));

[uu,vv] = meshgrid(param.us,param.vs);
[iu,iv] = meshgrid(1:param.nu,1:param.nv);

nmimg = numel(data3d);
szproj = [param.nu,param.nv];

data3d(isnan(data3d))=0;
data3d = permute(data3d,[1 3 2]);

[xx,zz] = meshgrid(param.xs,param.zs);
dist = sqrt((param.DSD)^2 + uu.^2 + vv.^2)./(param.DSD)*param.dy;
disttmp = dist';

M = sparse(prod(szproj),nmimg);

for iy = 1:param.ny

    Ratio = (param.ys(iy)+param.DSO)/(param.DSD);

    % if(Ratio<0)
    %     continue
    % end

    pu = uu*Ratio;
    pv = vv*Ratio;

    pu = (pu - xx(1,1))/(param.dx)+1;
    pv = (pv - zz(1,1))/(param.dz)+1;

    Indimg = sub2ind([param.nx,param.ny,param.nz],round(pu),iy*ones(size(pv)),round(pv));
    Indproj = sub2ind([param.nu,param.nv],iu,iv);
    M = M + sparse(Indproj,Indimg,disttmp(Indproj),prod(szproj),nmimg);

    tmp = (interp2((single(data3d(:,:,iy))),(single(pv)),(single(pu)),'nearest'));
    nnz(tmp)
    % figure(10);subplot(2,1,1);imshow(data3d(:,:,iy),[])
    % figure(10);subplot(2,1,2);imshow(tmp,[])

    tmp(isnan(tmp))=0;

    proj2d = proj2d + tmp';



end


proj2d = proj2d .* dist';

    test = permute(data3d,[1 3 2]);
    tmp2 = reshape(M*double(test(:)),[nu,nv]);
    figure(10);imshow([proj2d tmp2],[])
    pause(1)








