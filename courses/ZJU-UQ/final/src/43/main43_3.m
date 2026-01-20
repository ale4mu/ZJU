% 用于验证确定性问题中频率不影响峰值
clear; clc; close all;

% Example 4.3
cfg.gamma = 5/3;       % 绝热指数 5/3
cfg.L = 5.0;           % 区域长度 [0, 5]
cfg.T = 4.0;           % 终止时间
cfg.CFL = 0.6;
cfg.Nx = 5000;          % gPC-SG 的网格数
cfg.M = 5;             % gPC 阶数

fprintf('w = 1:');
num_points = 50;    
[nodes_sc, weights_sc] = lgwt(num_points, -1, 1);

Mean1 = zeros(1, cfg.Nx);
Var1  = zeros(1, cfg.Nx);
Mean2 = zeros(1, cfg.Nx);
Var2  = zeros(1, cfg.Nx);

dx = cfg.L / cfg.Nx;
x = (0.5 : cfg.Nx - 0.5) * dx;
parfor q = 1:num_points
    xi = nodes_sc(q);
    w_1 = 1 ; %  w(xi)
    w_2 = 1.1;
    [rho_final1, ~, ~] = collocation43(cfg, w_1);
    [rho_final2, ~, ~] = collocation43(cfg, w_2);

    w_q = weights_sc(q) * 0.5;
    Mean1 = Mean1 + rho_final1 * w_q;
    Var1  = Var1  + (rho_final1.^2) * w_q;
    Mean2 = Mean2 + rho_final2 * w_q;
    Var2  = Var2  + (rho_final2.^2) * w_q;
end
Std1 = sqrt(Var1 - Mean1.^2);
Std2 = sqrt(Var2 - Mean2.^2);


figure('Position', [100, 100, 1000, 400], 'Color', 'w');
subplot(1, 2, 1);

hold on;
plot(x, Mean1, 'k-', 'LineWidth', 1.5, 'DisplayName', 'w = 1');
hold off;
xlabel('x'); ylabel('Expectation');
title(['Expectation at t = ' num2str(cfg.T)]);
legend('Location', 'SouthEast');
grid on;
xlim([0, 5]);
ylim([0.99, 1.01]); 

subplot(1, 2, 2);
hold on;
plot(x, Mean2, 'k-', 'LineWidth', 1.5, 'DisplayName', 'w = 1.1');
hold off;
xlabel('x'); ylabel('Expectation');
title(['Expectation at t = ' num2str(cfg.T)]);
legend('Location', 'SouthEast');
grid on;
xlim([0, 5]);
ylim([0.99, 1.01]); 
