% 验证保双曲性
clear;clc;close all;

cfg.gamma = 1.4;
cfg.L = 1.0;
cfg.T = 0.2;
cfg.CFL = 0.6;
cfg.Nx = 100;
cfg.Ny = 100;
cfg.M = 4;


num_points = 20;
[nodes, weights] = lgwt(num_points, -1, 1);

Phi_val = zeros(cfg.M + 1, num_points);
for k = 0:cfg.M
    for q = 1:num_points
        Phi_val(k+1, q) = my_legendre(k, nodes(q));
    end
end
norm_sq = 2 ./ (2*(0:cfg.M)' + 1);

% 2D 网格
dx = cfg.L / cfg.Nx;
dy = cfg.L / cfg.Ny;
x_grid = (0.5:cfg.Nx-0.5) * dx;
y_grid = (0.5:cfg.Ny-0.5) * dy;
[X, Y] = meshgrid(x_grid, y_grid); % 注意 meshgrid 是 (Ny, Nx)
X = X';
Y = Y'; % 转置为 (Nx, Ny)

% 初始化 U_hat (4, M+1, Nx, Ny)
U_hat = zeros(4, cfg.M + 1, cfg.Nx, cfg.Ny);
for i = 1:cfg.Nx
    for j = 1:cfg.Ny
        x_val = x_grid(i);
        y_val = y_grid(j);
        % 计算该点的 gPC 系数
        U_coeffs = Initial_gPC(x_val, y_val, cfg, nodes, weights, Phi_val, norm_sq);
        U_hat(:, :, i, j) = U_coeffs;
    end
end


t = 0;
flag = true;
while t < cfg.T
    rho = squeeze(U_hat(1, 1, :, :));
    u   = squeeze(U_hat(2, 1, :, :)) ./ rho;
    v   = squeeze(U_hat(3, 1, :, :)) ./ rho;
    p   = (cfg.gamma-1)*(squeeze(U_hat(4, 1, :, :)) - 0.5*rho.*(u.^2+v.^2));
    c   = sqrt(cfg.gamma*abs(p)./rho);

    max_vel = max(max(abs(u)+c, abs(v)+c), [], 'all');
    fprintf("max_vel = %.4f\n", max_vel);
    dt = cfg.CFL * min(dx, dy) / max_vel;
    if t + dt > cfg.T, dt = cfg.T - t; end

    sqrt_13 = sqrt(13);
    c1 = (2 * dt) / (5 - sqrt_13 + sqrt(2 * (1 + sqrt_13)));
    c2 = ((7 + sqrt_13 - sqrt(2 * (1 + sqrt_13))) / 12) * dt;
    tau1 = c1;
    tau2 = c2;
    tau3 = (tau1^2) / (tau2 - tau1);
    tau4 = dt - (tau1 + tau2 + tau3);

    if abs(tau1 + tau2 + tau3 + tau4 - dt) > 1e-10
        warning('Time splitting coefficients do not sum to dt!');
    end
    
    U_hat = f_direction(U_hat,tau4, 'Y', dx, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    U_hat = f_direction(U_hat,tau3+tau4, 'X', dx, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    U_hat = f_direction(U_hat,tau3, 'Y', dx, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    U_hat = f_direction(U_hat, tau2, 'X', dy, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    U_hat = f_direction(U_hat,tau1 + tau2, 'Y', dx, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    U_hat = f_direction(U_hat, tau1, 'X', dy, cfg.gamma, Phi_val, weights, norm_sq, nodes,false);
    if t > 0.1 && flag
        fprintf('Checking Hyperbolicity at t = %.4f ...\n', t);
        
        ix = round(cfg.Nx / 2);
        iy = round(cfg.Ny / 2);
        
        U_slice = squeeze(U_hat(:, :, ix, iy));
        
        fprintf('  Checking X-Direction...\n');
        hyper_2D(U_slice, cfg.gamma, Phi_val, weights, norm_sq);

        fprintf('  Checking Y-Direction...\n');
        U_slice_Y = U_slice;
        U_slice_Y(2) = U_slice(3); % u <-> v
        U_slice_Y(3) = U_slice(2);
        hyper_2D(U_slice_Y, cfg.gamma, Phi_val, weights, norm_sq);
        
        flag = false;
        disp('Check Done. Press any key to continue...');
        pause; 
    end
    t = t + dt;
    fprintf('t = %.2f\n', t);
end













