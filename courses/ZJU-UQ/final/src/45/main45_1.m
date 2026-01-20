% 绘制等高线图和对角线图
clear;clc;close all;

cfg.gamma = 1.4;
cfg.L = 1.0;
cfg.T = 0.2;
cfg.CFL = 0.6;
cfg.Nx = 250;
cfg.Ny = 250;
cfg.M = 4;

[Mean_rho, Std_rho, X, Y] = gpc_SG_2D(cfg);
cfg_ref = cfg;
cfg_ref.Nx = 400;
cfg_ref.Ny = 400;
[Mean_Ref, Std_Ref] = collocaction45(cfg_ref);
fprintf('yes');

% 保存数据
save('e:\UQ\final\src\45\gpc_results_new.mat', 'Mean_rho', 'Std_rho', 'X', 'Y');
save('e:\UQ\final\src\45\ref_results_new.mat', 'Mean_Ref', 'Std_Ref');

% figure;
% contourf(X, Y, Mean_rho, 30);
% colorbar;
% title('Mean Density');

% figure;
% contourf(X, Y, Std_rho, 30);
% colorbar;
% title('Std Deviation');

% 等高线图
figure('Position', [100, 100, 1200, 500], 'Color', 'w');
subplot(1, 2, 1);
% contour(X, Y, Z, Levels)
contour(X, Y, Mean_rho, 30, 'LineWidth', 1.2); 
axis equal; axis([0 1 0 1]);
xlabel('x'); ylabel('y');
colorbar;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
contour(X, Y, Std_rho, 30, 'LineWidth', 1.2);
axis equal; axis([0 1 0 1]);
xlabel('x'); ylabel('y');
colorbar;
set(gca, 'FontSize', 12);


dx_ref = cfg.L / cfg_ref.Nx;
dy_ref = cfg.L / cfg_ref.Ny;
x_grid_ref = (0.5 : cfg_ref.Nx - 0.5) * dx_ref;
y_grid_ref = (0.5 : cfg_ref.Ny - 0.5) * dy_ref;
[X_Ref, Y_Ref] = meshgrid(x_grid_ref, y_grid_ref);
X_Ref = X_Ref'; Y_Ref = Y_Ref';

figure('Position', [100, 100, 1200, 500], 'Color', 'w');
subplot(1, 2, 1);
contour(X_Ref, Y_Ref, Mean_Ref, 30, 'LineWidth', 1.0);
axis equal; axis([0 1 0 1]);
xlabel('x'); ylabel('y');
colorbar;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
contour(X_Ref, Y_Ref, Std_Ref, 30, 'LineWidth', 1.0);
axis equal; axis([0 1 0 1]);
xlabel('x'); ylabel('y');
colorbar;
set(gca, 'FontSize', 12);
c = parcluster('local');

% 对角线图
s_diag = linspace(0, 1, 500); 
x_diag = s_diag * cfg.L;      % 实际 x 坐标
y_diag = s_diag * cfg.L;      % 实际 y 坐标

Mean_Diag = interp2(Y, X, Mean_rho, y_diag, x_diag, 'linear'); 
Std_Diag  = interp2(Y, X, Std_rho,  y_diag, x_diag, 'linear');
Mean_Ref_Diag = interp2(Y_Ref, X_Ref, Mean_Ref, y_diag, x_diag, 'linear');
Std_Ref_Diag = interp2(Y_Ref, X_Ref, Std_Ref, y_diag, x_diag, 'linear');
figure('Position', [100, 100, 1000, 400], 'Color', 'w');
subplot(1, 2, 1);
plot(s_diag(1:2:end), Mean_Diag(1:2:end), 'r*', 'LineWidth', 1.5,'MarkerSize', 4); 
hold on;
plot(s_diag, Mean_Ref_Diag, 'k-', 'LineWidth', 1.5);
%plot(s_diag(1:5:end), Mean_Diag(1:5:end), 'r*-', 'MarkerSize', 6);
grid on;
xlim([0, 1]);
ylim([0.2,1.1])

subplot(1, 2, 2);
plot(s_diag(1:2:end), Std_Diag(1:2:end), 'r*', 'LineWidth', 1.5,'MarkerSize', 4);
hold on;
plot(s_diag, Std_Ref_Diag, 'k-', 'LineWidth', 1.5);
%plot(s_diag(1:20:end), Std_Diag(1:20:end), 'r*-', 'MarkerSize', 6);
grid on;
xlim([0, 1]);