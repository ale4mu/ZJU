%% 主程序：gPC-Collocation方法求解随机ODE问题
clear; close all; clc;

N_list = [8, 16, 32, 64];      % 时间离散度列表
M_list = [1, 2, 3, 4, 5, 6];   % 随机配置点个数列表
beta = 1;                       % 初始条件
T = 1;                          % 时间区间[0,T]


fprintf('\n=============== 测试1: 固定M=10，变化N ===============\n');
fprintf('%-10s %-20s %-15s %-20s %-15s\n', ...
    'N', 'mean_error', 'mean_order', 'var_error', 'var_order');

past_mean_err = 0;
past_var_err = 0;
past_N = 0;

for i = 1:length(N_list)
    N = N_list(i);
    [mean_err, var_err] = solve_ode(N, 20);  % 固定M=10
    
    % 计算收敛阶
    if i == 1
        mean_order = 0;
        var_order = 0;
    else
        mean_order = order(past_mean_err, mean_err, past_N, N);
        var_order = order(past_var_err, var_err, past_N, N);
    end
    
    % 显示结果
    fprintf('%-10d %-20.6e %-15.4f %-20.6e %-15.4f\n', ...
        N, mean_err, mean_order, var_err, var_order);
    
    % 更新历史值
    past_mean_err = mean_err;
    past_var_err = var_err;
    past_N = N;
end


fprintf('\n=============== 测试2: 固定N=128，变化M ===============\n');
fprintf('%-15s %-20s %-20s\n', '多项式阶数M', 'mean_error', 'var_error');

% 准备保存结果
results_M = zeros(length(M_list), 3);  % [M, mean_err, var_err]

for i = 1:length(M_list)
    M = M_list(i);
    [mean_err, var_err] = solve_ode(128, M);  % 固定N=128
    
    % 显示结果
    fprintf('%-15d %-20.6e %-20.6e\n', M, mean_err, var_err);
    
    % 保存结果
    results_M(i, :) = [M, mean_err, var_err];
end


figure('Position', [100, 100, 1200, 800]);

semilogy(M_list, results_M(:, 2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(M_list, results_M(:, 3), 'b-d', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('随机配置点个数 M');
ylabel('误差');
title('误差随M的变化 (N=128)');
legend('期望误差', '方差误差', 'Location', 'best');
grid on;

