%% 主程序：gPC-Collocation方法求解随机PDE问题
clear; close all; clc;

% 参数设置
N_list = [32, 64, 128];        % 空间离散网格数
M_list = [1, 2, 3, 4, 5, 6]; % 随机配置点个数
T = 5.0;                       % 模拟结束时间 (对应C++中的 5.0)

% 存储结果
results = zeros(length(M_list), 3); % [M, mean_error, var_error]

fprintf('\n=============== 测试1: 固定M=10，变化N ===============\n');
fprintf('%-10s %-20s %-15s %-20s %-15s\n', ...
    'N', 'mean_error', 'mean_order', 'var_error', 'var_order');
for i = 1:length(N_list)
    N = N_list(i);
    [mean_err, var_err] = solve_pde(N,10, T);  % 固定M=10
    
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

fprintf('=============== 开始计算 PDE (T=%.1f) ===============\n', T);
fprintf('%-5s %-15s %-15s\n', 'M', 'Mean Error', 'Var Error');

% 这里固定 N=128 来测试 M 的收敛性 (对应C++的逻辑)
fixed_N = 4096; 

for i = 1:length(M_list)
    M = M_list(i);
    [mean_err, var_err] = solve_pde(fixed_N, M, T);
    
    results(i, :) = [M, mean_err, var_err];
    fprintf('%-5d %-15.6e %-15.6e\n', M, mean_err, var_err);
end
% 绘图
figure;
semilogy(results(:,1), results(:,2), '-bo', 'LineWidth', 1.5, 'DisplayName', '期望误差');
hold on;
semilogy(results(:,1), results(:,3), '-rs', 'LineWidth', 1.5, 'DisplayName', '方差误差');
grid on;
xlabel('配置点个数 M');
ylabel('误差');
legend('Location', 'best');

