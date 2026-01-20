% 试图验证谱收敛但误差完全被时间和空间离散误差主导
clear;clc;close all;

cfg.gamma = 1.4;      
cfg.L = 1.0;           % 区域长度
cfg.T = 0.2;           % 最终时间
cfg.CFL = 0.6;        
cfg.Nx = 100;         
cfg.Ny = 100;          
M_values = [0,1,2,3,4];

ref_data_path = 'e:\UQ\final\src\45\ref_results.mat';
load(ref_data_path);  % 加载Mean_Ref, Std_Ref

dx_sg = cfg.L / 100; 
x_sg = (0.5 : 100 - 0.5) * dx_sg;
y_sg = (0.5 : 100 - 0.5) * dx_sg;
[X_SG, Y_SG] = meshgrid(x_sg, y_sg);
X_SG = X_SG'; Y_SG = Y_SG';   

cfg_ref = cfg;
cfg_ref.Nx = 600;
cfg_ref.Ny = 600;
dx_ref = cfg.L / cfg_ref.Nx;
dy_ref = cfg.L / cfg_ref.Ny;
x_grid_ref = (0.5 : cfg_ref.Nx - 0.5) * dx_ref;
y_grid_ref = (0.5 : cfg_ref.Ny - 0.5) * dy_ref;

Mean_Ref_Down = interp2(y_grid_ref, x_grid_ref, Mean_Ref, y_sg, x_sg', 'cubic');
Std_Ref_Down  = interp2(y_grid_ref, x_grid_ref, Std_Ref,  y_sg, x_sg', 'cubic');

L2_errors_mean = zeros(size(M_values));
L2_errors_std = zeros(size(M_values));

for idx = 1:length(M_values)
    M = M_values(idx);
    cfg.M = M;
    fprintf('正在计算M = %d的gPC-SG解...\n', M);
    
    [Mean_rho, Std_rho, ~, ~] = gpc_SG_2D(cfg);

    L2_errors_mean(idx) = L2Error(Mean_rho, Mean_Ref_Down);
    L2_errors_std(idx) = L2Error(Std_rho, Std_Ref_Down);
end

figure('Position', [100, 100, 1000, 500], 'Color', 'w');

subplot(1,2,1);
semilogy(M_values, L2_errors_mean, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
title('gPC-SG均值误差随M的收敛图');
xlabel('gPC阶数M');
ylabel('log(L2误差)');
grid on;
set(gca, 'FontSize', 12);

subplot(1,2,2);
semilogy(M_values, L2_errors_std, 's--', 'LineWidth', 1.5, 'MarkerSize', 8);
title('gPC-SG标准差误差随M的收敛图');
xlabel('gPC阶数M');
ylabel('log(L2误差)');
grid on;
set(gca, 'FontSize', 12);
