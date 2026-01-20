function [Mean_rho,Std_rho] = gpc_SG4_1(cfg)

fprintf('M = %d \n',cfg.M);
num_points = 100;
[nodes,weights] = lgwt(num_points, -1, 1); % Legendre-Gauss 积分点和权重

% Phi_val: (M+1,num_points)
Phi_val = zeros(cfg.M + 1, num_points);
for k = 0:cfg.M
    for q = 1:num_points
        Phi_val(k+1, q) = my_legendre(k, nodes(q)); % k阶
    end
end

norm_sq = 2 ./ (2*(0:cfg.M)' + 1);

dx = cfg.L / cfg.Nx;
x = (0.5:cfg.Nx-0.5) * dx; % 单元中心
U_hat = zeros(3, cfg.M + 1, cfg.Nx);

% 初始化
for i = 1:cfg.Nx
    x_val = x(i);
    % 计算gPC系数    
    for q = 1:num_points
        xi_val = nodes(q);
        w = weights(q);

        rho_val = 1 + 0.2 * sin(2 * pi * x_val);
        u_val   = 0.8 + 0.2 * xi_val;
        p_val   = 1.0;
        E_val   = p_val / (cfg.gamma - 1) + 0.5 * rho_val * u_val^2;
        
        % U:[rho, rho*u, E]
        U_vec = [rho_val; rho_val * u_val; E_val];

        for k = 0:cfg.M
            basis = Phi_val(k+1, q);
            U_hat(:, k+1, i) = U_hat(:, k+1, i) + w * U_vec * basis;
        end
    end

    for k = 0:cfg.M
        U_hat(:, k+1, i) = U_hat(:, k+1, i) / norm_sq(k+1);
    end
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
    dt = 0.3 * dx^(5/3); % 验证5阶收敛阶用这个(main41_2.m)
    if t + dt > cfg.T, dt = cfg.T - t; end
   
    U_hat = RK3_SG(U_hat, dt, dx, cfg.gamma, Phi_val,weights,norm_sq,t,nodes);
    
    t = t + dt;
    if mod(t, 0.05) < dt, fprintf('t = %.4f\n', t); end
end

Mean_rho = squeeze(U_hat(1, 1, :))'; % 第0个系数

Var_rho = zeros(1, cfg.Nx);
for k = 1:cfg.M
    % k对应索引 k+1
    u_k = squeeze(U_hat(1, k+1, :))';
    norm = 1 / (2*k + 1); % 归一化
    Var_rho = Var_rho + u_k.^2 * norm;
end
Std_rho = sqrt(Var_rho);

end
