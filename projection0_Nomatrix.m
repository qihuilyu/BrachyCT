function proj2d = projection0_Nomatrix(param, img)

angle_rad = 0;
proj2d = (zeros(param.nu,param.nv,'single'));

[uu,vv] = meshgrid(param.us,param.vs);
[xx,yy] = meshgrid(param.xs,param.ys);

rx = (((xx.*cos(angle_rad) - yy.*sin(angle_rad)) - xx(1,1))/param.dx + 1);
ry = (((xx.*sin(angle_rad) + yy.*cos(angle_rad)) - yy(1,1))/param.dy + 1);


for iz = 1:param.nz

    img(:,:,iz) = interp2(img(:,:,iz),rx,ry, 'linear');

end

img(isnan(img))=0;
img = permute(img,[1 3 2]);

[xx,zz] = meshgrid(param.xs,param.zs);


for iy = 1:param.ny

    Ratio = (param.ys(iy)+param.DSO)/(param.DSD);

    if(Ratio<0)
        continue
    end

    pu = uu*Ratio;
    pv = vv*Ratio;

    pu = (pu - xx(1,1))/(param.dx)+1;
    pv = (pv - zz(1,1))/(param.dz)+1;

    tmp = (interp2((single(img(:,:,iy))),(single(pv)),(single(pu)),'linear'));

    % figure(10);subplot(2,1,1);imshow(data3d(:,:,iy),[])
    % figure(10);subplot(2,1,2);imshow(tmp,[])

    tmp(isnan(tmp))=0;

    proj2d = proj2d + tmp';
end

dist = sqrt((param.DSD)^2 + uu.^2 + vv.^2)./(param.DSD)*param.dy;

proj2d = proj2d .* dist';





