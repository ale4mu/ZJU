function [Mean_rho, Std_rho, x_grid] = gpc_SG4_3(cfg)
    % Example 4.3
    num_points = 150; % 积分点数
    [nodes, weights] = lgwt(num_points, -1, 1);
    
    Phi_val = zeros(cfg.M + 1, num_points);
    for k = 0:cfg.M
        for q = 1:num_points
            Phi_val(k+1, q) = my_legendre(k, nodes(q)); 
        end
    end

    norm_sq = 2 ./ (2*(0:cfg.M)' + 1);
    dx = cfg.L / cfg.Nx;
    x_grid = (0.5:cfg.Nx-0.5) * dx;
    
    % 初始条件: rho=1, u=0, p=0.6
    U_hat = zeros(3, cfg.M + 1, cfg.Nx);
    % rho_0 = 1
    U_hat(1, 1, :) = 1.0; 
    % m_0 = rho*u = 0 (已默认为0)
    
    p_init = 0.6;
    E_init = p_init / (cfg.gamma - 1);
    U_hat(3, 1, :) = E_init; 
    
    t = 0;
    while t < cfg.T
       
        rho0 = squeeze(U_hat(1, 1, :))';
        m0   = squeeze(U_hat(2, 1, :))';
        E0   = squeeze(U_hat(3, 1, :))';
        
        u0 = m0 ./ rho0;
        p0 = (cfg.gamma - 1) * (E0 - 0.5 * rho0 .* u0.^2);
        c0 = sqrt(cfg.gamma * abs(p0) ./ rho0);
        max_vel = max(abs(u0) + c0);
        
        dt = cfg.CFL * dx/max_vel;
        fprintf('dt = %.4f\n', dt);

        if t + dt > cfg.T, dt = cfg.T - t; end %走完最后一个时间步
        U_hat = RK3_SG(U_hat, dt, dx, cfg.gamma, Phi_val,weights, norm_sq, t,nodes);
        t = t + dt;
        if mod(t, 0.5) < dt, fprintf('t = %.4f / %.1f\n', t, cfg.T); end
    end
    Mean_rho = squeeze(U_hat(1, 1, :))'; 

    Var_rho = zeros(1, cfg.Nx);
    for i = 1:cfg.Nx
        coeffs = U_hat(1, :, i); % (M+1)
        u_vals = coeffs * Phi_val; % (1, N_quad)
        mean_val = Mean_rho(i);
        var_val = sum( (u_vals - mean_val).^2 .* weights' ) * 0.5;
        Var_rho(i) = var_val;
    end
    Std_rho = sqrt(Var_rho);
end
