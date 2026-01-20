% Example 4.1
function res = solve_SC(cfg,initial_func_handle)
%   cfg: PDE的参数，包含gamma, L, Nx, T, CFL
%   initial_func_handle: 初始化U的函数
    dx = cfg.L / cfg.Nx;
    x = (0.5:cfg.Nx-0.5) * dx; % 网格中心坐标
    %U = initial_func_handle(x);

    U = zeros(3, cfg.Nx);
    [g_nodes, g_weights] = lgwt(3, -1, 1);
    for i = 1:cfg.Nx
        x_center = (i-0.5)*dx;
        for g = 1:3
            x_val = x_center + 0.5*dx*g_nodes(g);
            w = g_weights(g) * 0.5;
            U_val = initial_func_handle(x_val); 
            U(:, i) = U(:, i) + w * U_val;
        end
    end

    t = 0;
    while t < cfg.T
        rho = U(1,:);
        u = U(2,:) ./ rho;
        p = (cfg.gamma - 1) * (U(3,:) - 0.5 * rho .* u.^2);
        c = sqrt(cfg.gamma * p ./ rho);
        max_lamda = max(abs(u) + c);

        dt = cfg.CFL * dx / max_lamda;
        %dt = 0.3 * (dx^(5/3));
        if t + dt > cfg.T, dt = cfg.T - t; end % 最后一个时间步

        U = RK3_SC(U, dt, dx, cfg.gamma);
        t = t + dt;
    end
    res = U(1,:);
end