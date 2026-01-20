function [Mean_Ref, Std_Ref] = collocation42(cfg)

    num_points_1d = 10; % 每个维度 10 个点，总共 100 个点
    [nodes_1d, weights_1d] = lgwt(num_points_1d, -1, 1);
    
    [XI1, XI2] = meshgrid(nodes_1d, nodes_1d);
    [W1, W2] = meshgrid(weights_1d, weights_1d);
    

    nodes_2d = [XI1(:), XI2(:)]'; % (2, N_quad)
    weights_2d = W1(:) .* W2(:);  % (N_quad, 1)
    N_quad = length(weights_2d);
    
    Results = zeros(N_quad, cfg.Nx);
    
    parfor q = 1:N_quad
        xi1 = nodes_2d(1, q);
        xi2 = nodes_2d(2, q);
        
        current_init_func = @(x) initial_func_Ex42(x, cfg, xi1, xi2);
        
        Results(q, :) = solve_SC(cfg, current_init_func);
    end
    
    Mean_Ref = zeros(1, cfg.Nx);
    Mean_Rho2 = zeros(1, cfg.Nx);
    
    pdf_norm = 0.25; 
    
    for q = 1:N_quad
        w_q = weights_2d(q) * pdf_norm;
        res = Results(q, :);
        
        Mean_Ref = Mean_Ref + res * w_q;
        Mean_Rho2 = Mean_Rho2 + (res.^2) * w_q;
    end
    
    Std_Ref = sqrt(Mean_Rho2 - Mean_Ref.^2);
end

function U = initial_func_Ex42(x, cfg, xi1, xi2)
    % (rho, u, p)(x, 0, xi) = (1.0 + 0.2(1 + 0.1 sin(pi*xi1/4)) sin(2pi*x), 0.8 + 0.2*xi2, 1)
    
    rho = 1.0 + 0.2 * (1 + 0.1 * sin(pi * xi1 / 4)) * sin(2 * pi * x);
    u   = (0.8 + 0.2 * xi2) * ones(size(x));
    p   = ones(size(x));
    
    E = p / (cfg.gamma - 1) + 0.5 * rho .* u.^2;
    
    U = [rho; rho.*u; E];
end