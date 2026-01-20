% 初始条件生成
function U = initial_func4_1(x,cfg,xi)
    rho = 1 + 0.2 * sin(2 * pi * x);
    u   = (0.8 + 0.2 * xi) * ones(size(x));
    p   = ones(size(x));
    E   = p / (cfg.gamma - 1) + 0.5 * rho .* u.^2;
    U = [rho; rho.*u; E];
end

