function M = projection0(param)

nu = param.nu;
nv = param.nv;

[uu,vv] = meshgrid(param.us,param.vs);

nmimg = prod([param.nx,param.ny,param.nz]);
szproj = [nu,nv];

[xx,zz] = meshgrid(param.xs,param.zs);
dist = sqrt((param.DSD)^2 + uu.^2 + vv.^2)./(param.DSD)*param.dy;
disttmp = dist';

M = sparse(prod(szproj),nmimg);

parfor iy = 1:param.ny

    Ratio = (param.ys(iy)+param.DSO)/(param.DSD);

    if(Ratio<=0)
        continue
    end

    pu = uu*Ratio;
    pv = vv*Ratio;

    pu = (pu - xx(1,1))/(param.dx)+1;
    pv = (pv - zz(1,1))/(param.dz)+1;
    [iu,iv] = meshgrid(1:param.nu,1:param.nv);

    ind_keep = find(~(round(pu)<1|round(pv)<1|round(pu)>param.nx|round(pv)>param.nz));
    pu = pu(ind_keep);
    pv = pv(ind_keep);
    iu = iu(ind_keep);
    iv = iv(ind_keep);

    Indimg = sub2ind([param.nx,param.ny,param.nz],round(pu),iy*ones(size(pv)),round(pv));
    Indproj = sub2ind([param.nu,param.nv],iu,iv);
    M = M + sparse(Indproj,Indimg,disttmp(Indproj),prod(szproj),nmimg);

end






