function [x_3d, Maincost] = IterRecon_TV_FISTA_BrachyCT (x0, A, g1,g2,g3,g4, gamma, mu, ImageSize)
% CT Iterative Reconstruction using BM3D Regularization and FISTA algorithm.
% Solves the optimization (1/2)(||A*W*x-g||_2)^2+gamma*||D*x||_1^{mu}.
% g are the Line Integral images. x is the 3D reconstructed image, or the optimization variable.
% AWx are a list of projection images, which are of the same dimension as g.
% W is a matrix to deal with boundary condition problems.
fprintf('\n\n\n*******Iterative Reconstruction*******\n\n');

% Initialization.
N = prod(ImageSize);
% xk = zeros(N,1); vk = xk;
xk = x0(:); vk = xk;
tk = 1e-06; % Initial step sizes.
r1 = 0.8; r2 = 0.1;

maxIter = 50; count = 1;
Maincost = zeros(maxIter,1);

% Iterative loop
fprintf('Begin Optimization...\n\n');
tic

while (count<=maxIter)
    t = tk/r1;
    subiteration = 0;
    if count == 2
        r2 = 0.5;
    end
    while 1 % Lipschitz Gradient Condition
        if(count == 1)
            theta = 1;
        else
            r = tk/t;
            theta = thetak*(sqrt(4*r+thetak^2) - thetak)/(2*r);
        end

        y = (1-theta)*xk + theta*vk;

        FPy = ForwardProj(y,A);
        DataDiff_y = FPy-[g1;g2;g3;g4];
        BPz1 = BackwardProj(DataDiff_y, A,ImageSize);
        DVP_y = DVP_Dx (mu, reshape(y,ImageSize));
        delta_fy = BPz1(:) + gamma/mu*(DVP_y(:));

        x = min(max(0,y-t*delta_fy),0.2);

        % Compute f(y)
        Dy = applyD3D (reshape(y,ImageSize));
        Huber_Dy_in = 1/(2*mu)*(Dy).^2;
        Huber_Dy_out = abs(Dy)-mu/2;
        Huber_Dy = Huber_Dy_in;
        Huber_Dy(abs(Dy)>mu) = Huber_Dy_out(abs(Dy)>mu);
        fy = 1/2*sum(DataDiff_y(:).^2) + gamma*sum(Huber_Dy(:));

        UpperBound_x = fy + sum(delta_fy.*(x-y)) + 1/(2*t)*sum((x-y).^2);

        % Compute f(x)
        FPx = ForwardProj(x, A);    % Para Version
        DataDiff_x = FPx-[g1;g2;g3;g4];

        Dx = applyD3D (reshape(x,ImageSize));
        Huber_Dx_in = 1/(2*mu)*(Dx).^2;
        Huber_Dx_out = abs(Dx)-mu/2;
        Huber_Dx = Huber_Dx_in;
        Huber_Dx(abs(Dx)>mu) = Huber_Dx_out(abs(Dx)>mu);

        fx = 1/2*sum(DataDiff_x(:).^2) + gamma*sum(Huber_Dx(:));

        if fx <= UpperBound_x
            break
        end

        t = r2*t
        subiteration = subiteration+1;
    end

    tk = t;
    thetak = theta;
    vk = xk + 1/theta*(x - xk);
    xk = x;

    Maincost(count) = fx;

    if mod(count,1) == 0
        opttime=toc
        figure(3);semilogy(Maincost,'r');hold on;title(['TV iteration:' num2str(count) ', gamma: ' num2str(gamma)]); pause(0.01);
        x_3d = reshape(x,ImageSize);
        figure(6);imshow(x_3d(:,:,ceil(end/2)), []);title(['TV iteration:' num2str(count) ', gamma: ' num2str(gamma)]); pause(0.01);
        figure(7);imshow(squeeze(x_3d(:,ceil(end/2),:)),[]);
    end
    if mod(count,50) == 0
        save(['TV_mu_' num2str(mu) '_gamma_' num2str(gamma) '_count_' num2str(count) '.mat'],'x_3d','Maincost','gamma','mu','count','tk','thetak','xk','vk','-v7.3');
    end
    count = count+1;
end

fprintf('\n\nEnd Optimization.');
fprintf('\n\n*******Iterative Reconstruction*******\n\n');

x_3d = reshape(x,ImageSize);

end



function Ax0 = ForwardProj(x0,A)

x1 = x0;
y1_all = A*x1(:);
x2 = permute(x0,[2,1,3]);
y2_all = A*x2(:);
x3 = flip(x0,2);
y3_all = A*x3(:);
x4 = flip(x2,2);
y4_all = A*x4(:);

Ax0 = [y1_all; y2_all; y3_all; y4_all];

end


function ATy = BackwardProj(y0,A,ImageSize)

y1 = y0(1:end/4);
y2 = y0(end/4+1:end/4*2);
y3 = y0(end/4*2+1:end/4*3);
y4 = y0(end/4*3+1:end);

ATy1 = reshape(A'*y1,ImageSize);
ATy2 = reshape(A'*y2,ImageSize);
ATy2 = permute(ATy2,[2,1,3]);
ATy3 = reshape(A'*y3,ImageSize);
ATy3 = flip(ATy3,2);
ATy4 = reshape(A'*y4,ImageSize);
ATy4 = permute(ATy4,[2,1,3]);
ATy4 = flip(ATy4,2);

ATy = ATy1+ATy2+ATy3+ATy4;

end




