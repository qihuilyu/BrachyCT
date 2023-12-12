function vol = backprojection(proj,param,iview)

angle_rad = param.deg(iview)/360*2*pi;
vol = zeros(param.nx,param.ny,param.nz,'single');

[xx,yy] = meshgrid(param.xs,param.ys);

rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);

pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));

[uu,vv] = meshgrid(param.us,param.vs);
dist = sqrt((param.DSD)^2 + uu.^2 + vv.^2)./(param.DSD)*param.dy;


for iz = 1:param.nz
    pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
    vol(:,:,iz) = (interp2(dist.*proj',pu,pv,param.interptype));
end

vol(isnan(vol))=0;

return

















