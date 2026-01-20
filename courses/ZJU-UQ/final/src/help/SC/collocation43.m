% Example 4.3的随机配点法
function [rho, u, p] = collocation43(cfg, w_freq)
   
    Nx = cfg.Nx;
    dx = cfg.L / Nx;
    gamma = cfg.gamma;
    
    rho = ones(1, Nx);
    u   = zeros(1, Nx);
    p   = 0.6 * ones(1, Nx);
    E   = p/(gamma-1) + 0.5*rho.*u.^2;
    
    U = [rho; rho.*u; E]; % (3, Nx)
    
    t = 0;
    while t < cfg.T
        % if any(U(1, :) < 1e-4) || any(U(3, :) < 0)
        %     fprintf('Warning: Fixing non-physical state at t=%f\n', t);
        
        %     mask_rho = U(1, :) < 1e-2;
        %     U(1, mask_rho) = 1e-2;
        %     U(2, mask_rho) = 0; % 动量清零，防止速度过大
            
        %     % 修正能量/压力
        %     rho_fixed = U(1, :);
        %     u_fixed = U(2, :) ./ rho_fixed;
        %     p = (gamma-1) * (U(3, :) - 0.5 * rho_fixed .* u_fixed.^2);
            
        %     mask_p = p < 1e-2;
        %     % 重建能量以保证最小压力
        %     p(mask_p) = 1e-2;
        %     U(3, :) = p/(gamma-1) + 0.5 * rho_fixed .* u_fixed.^2;
        % end
        rho = U(1,:); 
        u = U(2,:)./rho; 
        p = (gamma-1)*(U(3,:) - 0.5*rho.*u.^2);
        c = sqrt(gamma*p./rho);
        lambda_max = max(abs(u)+c);
        dt = cfg.CFL * dx / lambda_max;
        if t+dt > cfg.T, dt = cfg.T - t; end
        
        U0 = U;

        L1 = upwind_SC43(U0, t, w_freq, dx, gamma);
        U1 = U0 + dt * L1;
        
        L2 = upwind_SC43(U1, t+dt, w_freq, dx, gamma);
        U2 = 0.75*U0 + 0.25*(U1 + dt*L2);
        
        L3 = upwind_SC43(U2, t+0.5*dt, w_freq, dx, gamma);
        U = (1/3)*U0 + (2/3)*(U2 + dt*L3);
        
        t = t + dt;
    end
    
    rho = U(1,:);
    u = U(2,:)./rho;
    p = (gamma-1)*(U(3,:) - 0.5*rho.*u.^2);
end