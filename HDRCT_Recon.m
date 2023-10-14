clear
close all


patname = 'T2';
datatestpath = fullfile('/Users/qihuilyu/Desktop/Data/HDRCT/Datatest/',patname);
load(fullfile(datatestpath,'Datatest.mat'));
imgdx = spatial.PixelSpacings(1);
temp = diff(spatial.PatientPositions);
imgdz = temp(1,3);
imgsz = size(CTimage);

imgdx = spatial.PixelSpacings(1);
temp = diff(spatial.PatientPositions);
imgdz = temp(1,3);
[CenterOfMass, FOV, deltax, deltay, deltaz]  = GetDwellPos(Mask_prostate,4,4,5,imgdx,imgdx,imgdz);

count = 1;
for dz = deltaz
    for dx = deltax
        for dy = deltay
            cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, 'dfs', inf, ...
                'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
                'nt', 1024, 'dt', 500/1024, ...
                'na', 1, 'orbit', 360,'orbit_start',0,...
                'pitch',0, 'source_z0', dz);

            ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
                'dx', imgdx, 'dz', imgdz, ...
                'offset_x', dy, 'offset_y', 0, 'offset_z', 0);

%             figure; cg.plot([ig])

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
            %         figure(1);imshow(y_cpu',[])
            %         pause(2)
            count = count+1;
        end
    end
end



figure;gif('Projections.gif')
for ii = 1:numel(y_cpus)
    y_all(:,:,ii) = y_cpus{ii};
    imshow(y_cpus{ii},[12,30]);
    gif
end




for dz = deltaz
    for dx = deltax
        for dy = deltay
            cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, 'dfs', inf, ...
                'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
                'nt', 1024, 'dt', 500/1024, ...
                'na', 1, 'orbit', 360,'orbit_start',90,...
                'pitch',0, 'source_z0', dz);

            ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
                'dx', imgdx, 'dz', imgdz, ...
                'offset_x', dy, 'offset_y', 0, 'offset_z', 0);

%             figure; cg.plot([ig])

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
            %         figure(1);imshow(y_cpu',[])
            %         pause(2)
            count = count+1;
        end
    end
end

for dz = deltaz
    for dx = deltax
        for dy = deltay
            cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, 'dfs', inf, ...
                'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
                'nt', 1024, 'dt', 500/1024, ...
                'na', 1, 'orbit', 360,'orbit_start',180,...
                'pitch',0, 'source_z0', dz);

            ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
                'dx', imgdx, 'dz', imgdz, ...
                'offset_x', dy, 'offset_y', 0, 'offset_z', 0);

%             figure; cg.plot([ig])
            dx
            dy
            dz
            x0 = flip(flip(LAC_380keV,1),2);
            %         x0(204:207, 51:54, 86) = 10;
            %          figure;imshow3D(x0,[])
            % figure;imshow(flip(x0(:,:,100),1),[])
            f.ptype = 'sf2';
            A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
            A_cpu.arg.mexfun = @jf_mex;
            A_cpus{count} = A_cpu;
            y_cpu = A_cpu*x0;
            y_cpus{count} = y_cpu;
            %         figure(1);imshow(y_cpu',[])
            %         pause(2)

            count = count+1;
        end
    end
end

for dz = deltaz
    for dx = deltax
        for dy = deltay
            cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, 'dfs', inf, ...
                'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
                'nt', 1024, 'dt', 500/1024, ...
                'na', 1, 'orbit', 360,'orbit_start',270,...
                'pitch',0, 'source_z0', dz);

            ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
                'dx', imgdx, 'dz', imgdz, ...
                'offset_x', dy, 'offset_y', 0, 'offset_z', 0);

%             figure; cg.plot([ig])

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
            %         figure(1);imshow(y_cpu',[])
            %         pause(2)
            count = count+1;
        end
    end
end


% gif('myfile.gif','DelayTime',1/24,'resolution',400)
%
%% Write the rest of the frames 
% After the first frame has been written, write each subsequent frame simply by calling |gif| without
% any options. Here we loop through the remaining 29 frames:
% 
%  for k = 2:29
%     set(h,'Zdata',Z*t(k))
%     gif
%  end


% ForBack.applyFP = @(x) for ii = 1:numel(A_cpus) A{ii}*x end;
% ForBack.applyBP = @(x) A'*x;
% 
% 
% gamma = 500;
% mu = 1e-05;
% 
% [x_TV, Maincost_TV] = IterRecon_TV_FISTA (y_all(:), gamma, mu, [ig.nx, ig.ny, ig.nz]);
% 


% x0 = zeros(ig.nx, ig.ny, ig.nz);
% x0(180:210,250:262,:) = 1;
% figure;imshow(flip(x0(:,:,100),1),[])
% f.ptype = 'sf2';
% A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
% A_cpu.arg.mexfun = @jf_mex;
% y_cpu = A_cpu*x0;
% figure;imshow(y_cpu,[])


