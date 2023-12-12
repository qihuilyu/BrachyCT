function [x_3d, Maincost] = IterRecon_TV_FISTA_Simple_2D (Ind, A, g, gamma, mu, ImageSize)
% CT Iterative Reconstruction using BM3D Regularization and FISTA algorithm.
% Solves the optimization (1/2)(||A*W*x-g||_2)^2+gamma*||D*x||_1^{mu}.
% g are the Line Integral images. x is the 3D reconstructed image, or the optimization variable.
% AWx are a list of projection images, which are of the same dimension as g.
% W is a matrix to deal with boundary condition problems.
fprintf('\n\n\n*******Iterative Reconstruction*******\n\n');

g0 = g;
g0(Ind,:) = 0;
g = g(:);


% Initialize the input image
g(g<0)=0;
g(logical(isnan(g)))=0;

% Initialization.
N = prod(ImageSize);
xk = zeros(N,1); vk = xk;
% xk = x0(:); vk = xk;
tk = 1e-06; % Initial step sizes.
r1 = 0.8; r2 = 0.1;

maxIter = 200; count = 1;
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
        
        FPy = A*y;
        FPy = reshape(FPy,size(g0));
        FPy(Ind,:) = 0;
        DataDiff_y = FPy-g0;
        BPz1 = A'*DataDiff_y(:); 
        DVP_y = DVP_Dx_2D (mu, reshape(y,ImageSize));
        delta_fy = BPz1(:) + gamma/mu*(DVP_y(:));
        
        x = min(max(0,y-t*delta_fy),0.2);
        
        % Compute f(y)
        Dy = applyD2D (reshape(y,ImageSize));
        Huber_Dy_in = 1/(2*mu)*(Dy).^2;
        Huber_Dy_out = abs(Dy)-mu/2;
        Huber_Dy = Huber_Dy_in;
        Huber_Dy(abs(Dy)>mu) = Huber_Dy_out(abs(Dy)>mu);
        fy = 1/2*sum(DataDiff_y(:).^2) + gamma*sum(Huber_Dy(:));
        
        UpperBound_x = fy + sum(delta_fy.*(x-y)) + 1/(2*t)*sum((x-y).^2);
        
        % Compute f(x)
        FPx = A*x;    % Para Version
        DataDiff_x = FPx-g;
        DataDiff_x = reshape(DataDiff_x,size(g0));
        DataDiff_x(Ind,:) = 0;
        
        Dx = applyD2D (reshape(x,ImageSize));
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
    
    if mod(count,5) == 0
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


