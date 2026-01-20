% 用于SG方法计算的rho的期望和标准差与参考解的对比
% 参考解使用随机配点法+1阶迎风格式计算
clear; clc; close all;

% Example 4.3
cfg.gamma = 5/3;       % 绝热指数 5/3
cfg.L = 5.0;           % 区域长度 [0, 5]
cfg.T = 4.0;           % 终止时间
cfg.CFL = 0.6;
cfg.Nx = 100;          % gPC-SG 的网格数
cfg.M = 5;             % gPC 阶数

fprintf('Running gPC-SG Method (M=%d, Nx=%d)...\n', cfg.M, cfg.Nx);
tic;
[Mean_SG, Std_SG, x_SG] = gpc_SG4_3(cfg);
t_sg = toc;
fprintf('gPC-SG Done. Time: %.2f s\n', t_sg);



fprintf('Running Reference Solution (SC, Nx=1000)...\n');
cfg_ref = cfg;
cfg_ref.Nx = 5000;     
num_sc_points = 50;    % 配点数
[nodes_sc, weights_sc] = lgwt(num_sc_points, -1, 1);

Mean_Ref = zeros(1, cfg_ref.Nx);
Var_Ref  = zeros(1, cfg_ref.Nx);

dx_ref = cfg_ref.L / cfg_ref.Nx;
x_Ref = (0.5 : cfg_ref.Nx - 0.5) * dx_ref;
tic;
parfor q = 1:num_sc_points
    xi = nodes_sc(q);
    w_freq = 1 + 0.1 * xi; %  w(xi)
    [rho_final, ~, ~] = collocation43(cfg_ref, w_freq);

    w_q = weights_sc(q) * 0.5;
    Mean_Ref = Mean_Ref + rho_final * w_q;
    Var_Ref  = Var_Ref  + (rho_final.^2) * w_q;
end
Std_Ref = sqrt(Var_Ref - Mean_Ref.^2);
fprintf('Mean Density: %.4f\n', Mean_Ref(1,10));
t_ref = toc;
fprintf('Reference Solution Done. Time: %.2f s\n', t_ref);


figure('Position', [100, 100, 1000, 400], 'Color', 'w');
subplot(1, 2, 1);

hold on;
plot(x_Ref, Mean_Ref, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Reference (SC)');
plot(x_SG, Mean_SG, 'r*', 'MarkerSize', 4, 'LineWidth', 1.0, 'DisplayName', 'gPC-SG (M=5)');
hold off;
xlabel('x'); ylabel('Expectation');
title(['Expectation at t = ' num2str(cfg.T)]);
grid on;
xlim([1, 5]);
ylim([0.99, 1.01]); 


subplot(1, 2, 2);
hold on;
plot(x_Ref, Std_Ref, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Reference (SC)');
plot(x_SG, Std_SG, 'r*', 'MarkerSize', 4, 'LineWidth', 1.0, 'DisplayName', 'gPC-SG (M=5)');
hold off;
xlabel('x'); ylabel('Std');
title(['Std at t = ' num2str(cfg.T)]);
legend('Location', 'NorthEast');
grid on;
xlim([1, 5]);
ylim([0, 0.007]); 

% sgtitle('Comparison of gPC-SG and Reference Solution (Example 4.3)');