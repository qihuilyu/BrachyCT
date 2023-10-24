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



[CenterOfMass, FOV, deltax, deltay, deltaz]  = GetDwellPos(Mask_prostate,4,4,5,imgdx,imgdx,imgdz);

dx = 0; dy = 0;
count = 1;




%%
for dz = deltaz
    % for dx = deltax
        for dy = deltay
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
            
            ATys{count} = A_cpu'*(y_cpu(:));

            %         figure(1);imshow(y_cpu',[])
            %         pause(2)
            count = count+1;
    %     end
    end
end

% figure;gif('Projections.gif')
% for ii = 1:numel(y_cpus)
%     y_all(:,:,ii) = y_cpus{ii};
%     imshow(y_cpus{ii},[12,30]);
%     gif
% end

% %%
% for dz = deltaz
%     for dx = deltax
%         for dy = deltay
%             cg = ct_geom('fan', 'dsd', 250+dx, 'dod', 250, 'dfs', inf, ...
%                 'ns', 1024, 'ds', 500/1024, 'offset_s', dy/(500/1024), ...
%                 'nt', 1024, 'dt', 500/1024, ...
%                 'na', 1, 'orbit', 360,'orbit_start',90,...
%                 'pitch',0, 'source_z0', dz);
% 
%             ig = image_geom('nx', imgsz(1), 'ny', imgsz(2), 'nz', imgsz(3), ...
%                 'dx', imgdx, 'dz', imgdz, ...
%                 'offset_x', dy, 'offset_y', 0, 'offset_z', 0);
% 
% %             figure; cg.plot([ig])
% 
%             x0 = LAC_380keV;
%             %         x0(204:207, 51:54, 86) = 10;
%             %          figure;imshow3D(x0,[])
%             % figure;imshow(flip(x0(:,:,100),1),[])
%             f.ptype = 'sf2';
%             A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
%             A_cpu.arg.mexfun = @jf_mex;
%             A_cpus{count} = A_cpu;
%             y_cpu = A_cpu*x0;
%             y_cpus{count} = y_cpu;
%             %         figure(1);imshow(y_cpu',[])
%             %         pause(2)
%             count = count+1;
%         end
%     end
% end
% 


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



% x0 = zeros(ig.nx, ig.ny, ig.nz);
% x0(180:210,250:262,:) = 1;
% figure;imshow(flip(x0(:,:,100),1),[])
% f.ptype = 'sf2';
% A_cpu = Gcone(cg, ig, 'type', f.ptype, 'nthread', 40);
% A_cpu.arg.mexfun = @jf_mex;
% y_cpu = A_cpu*x0;
% figure;imshow(y_cpu,[])


