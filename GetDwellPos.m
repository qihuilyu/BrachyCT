function [CenterOfMass, FOV, deltax, deltay, deltaz] = GetDwellPos(PTV,Nx,Ny,dNz,dx,dy,dz)

ImageX=sum(sum(PTV,3),2);
MassSum=sum(ImageX);
Xc=0;Yc=0;Zc=0;
for i = 1:length(ImageX)
    Xc = Xc + i*ImageX(i)/MassSum;
end
Xc = floor(Xc+0.5);

ImageY=sum(sum(PTV,3),1);
for i = 1:length(ImageY)
    Yc = Yc + i*ImageY(i)/MassSum;
end
Yc = floor(Yc+0.5);

ImageZ=sum(sum(PTV,2),1);
for i = 1:length(ImageZ)
    Zc = Zc + i*ImageZ(i)/MassSum;
end
Zc = floor(Zc+0.5);

CenterOfMass = [Xc, Yc, Zc];

ind = find(ImageX);
FOV(1) = (max(ind) - min(ind))*dx;
deltax = (min(ind):(max(ind) - min(ind))/Nx:max(ind)) - CenterOfMass(1);
deltax = (deltax(2:end) - (max(ind) - min(ind))/Nx/2)*dx;

ind = find(ImageY);
FOV(2) = (max(ind) - min(ind))*dy;
deltay = (min(ind):(max(ind) - min(ind))/Ny:max(ind)) - CenterOfMass(2);
deltay = (deltay(2:end) - (max(ind) - min(ind))/Ny/2)*dy;

ind = find(ImageZ);
FOV(3) = (max(ind) - min(ind))*dz;
deltaz = (min(ind)- CenterOfMass(3))*dz:dNz:(max(ind)- CenterOfMass(3))*dz;

end
