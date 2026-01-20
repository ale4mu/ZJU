function rho_final = solve_2D_46(cfg, xi_val)
    % Example 4.6 
    % 密度扰动: rho = rho0 * (1 + 0.1 * xi)
    
    cfg_det = cfg;
    cfg_det.M = 0;
    U_hat = zeros(4, 1, cfg.Nx, cfg.Ny);

    dx = cfg.L / cfg.Nx;
    dy = cfg.L / cfg.Ny;
    x_grid = (0.5:cfg.Nx-0.5) * dx;
    y_grid = (0.5:cfg.Ny-0.5) * dy;

    for i = 1:cfg.Nx
        for j = 1:cfg.Ny
            x = x_grid(i); y = y_grid(j);

            if x > 0.5 && y > 0.5
                rho0=0.5197; u0=0.1; v0=0.1; p0=0.4;
            elseif x < 0.5 && y > 0.5
                rho0=1; u0=-0.6259; v0=0.1; p0=1;
            elseif x < 0.5 && y < 0.5
                rho0=0.8; u0=0.1; v0=0.1; p0=1;
            else
                rho0=1; u0=0.1; v0=-0.6259; p0=1;
            end
            
            rho = rho0 * (1 + 0.1 * xi_val);
            u = u0;
            v = v0;
            p = p0;

            U_hat(1, 1, i, j) = rho;
            U_hat(2, 1, i, j) = rho*u;
            U_hat(3, 1, i, j) = rho*v;
            E = p/(cfg.gamma-1) + 0.5*rho*(u^2+v^2);
            U_hat(4, 1, i, j) = E;
        end
    end

    nodes_dummy = 0;
    weights_dummy = 2;
    Phi_dummy = 1;
    norm_dummy = 2;

    t = 0;
    while t < cfg.T
       
        rho = squeeze(U_hat(1, 1, :, :));
        u   = squeeze(U_hat(2, 1, :, :)) ./ rho;
        v   = squeeze(U_hat(3, 1, :, :)) ./ rho;
        p   = (cfg.gamma-1)*(squeeze(U_hat(4, 1, :, :)) - 0.5*rho.*(u.^2+v.^2));
        c   = sqrt(cfg.gamma*abs(p)./rho);

        max_vel = max(max(abs(u)+c, abs(v)+c), [], 'all');
        dt = cfg.CFL * min(dx, dy) / max_vel;
        if t + dt > cfg.T, dt = cfg.T - t; end
        
        U_hat = f_direction(U_hat, 0.5*dt, 'X', dy, cfg.gamma, Phi_dummy, weights_dummy, norm_dummy, nodes_dummy, true);
        U_hat = f_direction(U_hat, dt,     'Y', dx, cfg.gamma, Phi_dummy, weights_dummy, norm_dummy, nodes_dummy, true);
        U_hat = f_direction(U_hat, 0.5*dt, 'X', dy, cfg.gamma, Phi_dummy, weights_dummy, norm_dummy, nodes_dummy, true);
        
        t = t + dt;
    end

    rho_final = squeeze(U_hat(1, 1, :, :));
end