function [Mean_rho,Std_rho] = gPC_SG4_2(cfg)

fprintf('M = %d \n',cfg.M);
num_points_1d = 20; 
[Phi_val, nodes_2d, weights_2d, norm_sq] = build_2D_gPC(cfg.M, num_points_1d);

[nModes, N_quad] = size(Phi_val); 

xi1 = nodes_2d(1, :)'; % (N_quad, 1)
xi2 = nodes_2d(2, :)'; % (N_quad, 1)
U_hat = zeros(3, nModes, cfg.Nx);

dx = cfg.L / cfg.Nx;
x = (0.5:cfg.Nx-0.5) * dx;

[g_nodes, g_weights] = lgwt(3, -1, 1);
for i = 1:cfg.Nx
    %x_val = x(i);
    x_center = (i-0.5) * dx;
    U_hat_cell = zeros(3, nModes);
    
    for g = 1:3
        % 映射到物理坐标
        x_val = x_center + 0.5 * dx * g_nodes(g);
        w_g = g_weights(g) * 0.5;

        % rho(xi)
        rho_phys = 1.0 + 0.2 * (1 + 0.1 * sin(pi * xi1 / 4)) * sin(2 * pi * x_val);
        
        % u(xi)
        u_phys = 0.8 + 0.2 * xi2;
        
        % p(xi)
        p_phys = ones(N_quad, 1);
        
        % E(xi)
        E_phys = p_phys / (cfg.gamma - 1) + 0.5 * rho_phys .* (u_phys.^2);
        
        % 组装 U_phys (3, N_quad)
        U_phys = [rho_phys'; (rho_phys .* u_phys)'; E_phys']; 
        
        for k = 1:nModes
            basis = Phi_val(k, :)'; % (N_quad, 1)
            % 动量 m
            %val_m = sum(U_phys .* weights_2d' .* basis', 2); % (3, 1)
            %U_hat(:, k, i) = val_m / norm_sq(k);
            val_m = sum(U_phys .* weights_2d' .* basis', 2);
            U_hat_cell(:, k) = U_hat_cell(:, k) + w_g * (val_m / norm_sq(k));
        end
    end
    U_hat(:, :, i) = U_hat_cell;
end

t = 0;
while t < cfg.T

    rho0 = squeeze(U_hat(1, 1, :))';
    m0   = squeeze(U_hat(2, 1, :))';
    E0   = squeeze(U_hat(3, 1, :))';

    u0   = m0 ./ rho0;
    p0   = (cfg.gamma - 1) * (E0 - 0.5 * rho0 .* u0.^2);
    c0   = sqrt(cfg.gamma * abs(p0) ./ rho0);
    max_vel = max(abs(u0) + c0);

    %dt = cfg.CFL * dx; % 验证谱收敛用这个(main41_1.m)
    dt = 0.2 * dx^(5/3); % 验证5阶收敛阶用这个(main41_2.m)
    if t + dt > cfg.T, dt = cfg.T - t; end

    U_hat = RK3_SG(U_hat, dt, dx, cfg.gamma, Phi_val,weights_2d,norm_sq,t,nodes_2d);

    t = t + dt;
    if mod(t, 0.05) < dt, fprintf('t = %.4f\n', t); end
end

Mean_rho = squeeze(U_hat(1, 1, :))'; % 第0个系数

Var_rho = zeros(1, cfg.Nx);
for k = 2:nModes 
    u_k = squeeze(U_hat(1, k, :))';
    Var_rho = Var_rho + u_k.^2 * norm_sq(k) * 0.25;
end
Std_rho = sqrt(Var_rho);

end
