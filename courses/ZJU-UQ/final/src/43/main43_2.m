% 用于验证保双曲性
clear; clc; close all;
cfg.gamma = 5/3;       % 绝热指数 5/3
cfg.L = 5.0;           % 区域长度 [0, 5]
cfg.T = 4.0;           % 终止时间
cfg.CFL = 0.6;
cfg.Nx = 100;          % gPC-SG 的网格数
cfg.M = 5;             % gPC 阶数
% Example 4.3
fprintf('Running Example 4.3 with M = %d, Nx = %d\n', cfg.M, cfg.Nx);

num_points = 150; % guass积分点数
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

p_init = 0.6;
E_init = p_init / (cfg.gamma - 1);
U_hat(3, 1, :) = E_init;

t = 0;
flag = true;
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

    if t + dt > cfg.T
        dt = cfg.T - t; 
        if flag
            fprintf('Verifying Hyperbolicity at t = %.4f ...\n', t);
            mid_idx = round(cfg.Nx / 2);
            U_hat_slice = U_hat(:, :, mid_idx);
            hyper(U_hat_slice, cfg.gamma, Phi_val, weights, norm_sq);
            flag = false;
            break;
        end
    end
    U_hat = RK3_SG(U_hat, dt, dx, cfg.gamma, Phi_val,weights, norm_sq, t,nodes);
    t = t + dt;
    if mod(t, 0.5) < dt, fprintf('t = %.4f / %.1f\n', t, cfg.T); end
end

