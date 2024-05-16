function M = projection0(param, img)

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

    % tmp = (interp2((squeeze(single(img(:,iy,:)))),(single(pv)),(single(pu)),'nearest'));
    % tmp2 = (interp2((squeeze(single(img(:,iy,:)))),(single(pv)),(single(pu)),'linear'));
    % figure(1);imshow([tmp tmp2],[])

    puind1 = floor(pu);
    puind2 = ceil(pu);
    puweight1 = puind2-pu;
    puweight2 = pu-puind1;
    puweight1(puind1==puind2) = 0.5;
    puweight2(puind1==puind2) = 0.5;

    pvind1 = floor(pv);
    pvind2 = ceil(pv);
    pvweight1 = pvind2-pv;
    pvweight2 = pv-pvind1;
    pvweight1(pvind1==pvind2) = 0.5;
    pvweight2(pvind1==pvind2) = 0.5;

    ind_keep = find(~(round(puind1)<1|round(pvind1)<1|round(puind2)>param.nx|round(pvind2)>param.nz));
    pu = pu(ind_keep);
    pv = pv(ind_keep);
    iu = iu(ind_keep);
    iv = iv(ind_keep);
    puind1 = puind1(ind_keep);
    puind2 = puind2(ind_keep);
    puweight1 = puweight1(ind_keep);
    puweight2 = puweight2(ind_keep);

    pvind1 = pvind1(ind_keep);
    pvind2 = pvind2(ind_keep);
    pvweight1 = pvweight1(ind_keep);
    pvweight2 = pvweight2(ind_keep);

    Indimg1 = sub2ind([param.nx,param.ny,param.nz],puind1,iy*ones(size(pv)),pvind1);
    Indimg2 = sub2ind([param.nx,param.ny,param.nz],puind1,iy*ones(size(pv)),pvind2);
    Indimg3 = sub2ind([param.nx,param.ny,param.nz],puind2,iy*ones(size(pv)),pvind1);
    Indimg4 = sub2ind([param.nx,param.ny,param.nz],puind2,iy*ones(size(pv)),pvind2);
    pweight1 = puweight1.*pvweight1;
    pweight2 = puweight1.*pvweight2;
    pweight3 = puweight2.*pvweight1;
    pweight4 = puweight2.*pvweight2;

    Indproj = sub2ind([param.nu,param.nv],iu,iv);
    tmpM1 = sparse(Indproj,Indimg1,disttmp(Indproj).*pweight1,prod(szproj),nmimg);
    tmpM2 = sparse(Indproj,Indimg2,disttmp(Indproj).*pweight2,prod(szproj),nmimg);
    tmpM3 = sparse(Indproj,Indimg3,disttmp(Indproj).*pweight3,prod(szproj),nmimg);
    tmpM4 = sparse(Indproj,Indimg4,disttmp(Indproj).*pweight4,prod(szproj),nmimg);
    tmpMv2 = tmpM1+tmpM2+tmpM3+tmpM4;
    M = M + tmpMv2;

    % Indimg = sub2ind([param.nx,param.ny,param.nz],round(pu),iy*ones(size(pv)),round(pv));
    % Indproj = sub2ind([param.nu,param.nv],iu,iv);
    % tmpM = sparse(Indproj,Indimg,disttmp(Indproj),prod(szproj),nmimg);
    % M = M + tmpM;

    % tmp3 = reshape(tmpM*img(:),[param.nu,param.nv]);
    % tmp4 = reshape(tmpMv2*img(:),[param.nu,param.nv]);
    % figure(1);imshow([[tmp.*disttmp tmp2.*disttmp];[tmp3' tmp4'];[tmp.*disttmp tmp2.*disttmp]-[tmp3' tmp4']],[])
    % 

end






