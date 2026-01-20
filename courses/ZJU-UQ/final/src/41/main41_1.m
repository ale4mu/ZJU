% 用于验证谱收敛
clear;clc;close all;

cfg.gamma = 1.4;       % 气体常数（绝热率）
cfg.L = 1.0;           % 区域长度
cfg.T = 0.2;
cfg.CFL = 0.6;
cfg.Nx = 320; 

M_values = [0,1, 2, 3, 4, 5]; 
L2_errors_M = zeros(size(M_values));
Linf_errors_M = zeros(size(M_values));
L2_errors_std = zeros(size(M_values));
Linf_errors_std = zeros(size(M_values));

initial_func = @initial_func4_1;

for idx_M = 1:length(M_values)
    M = M_values(idx_M);
    cfg.M = M;
    
    dx = cfg.L / cfg.Nx;
    x_grid = (0.5:cfg.Nx-0.5) * dx;

    theta = 2 * pi * (x_grid - 0.8 * cfg.T);
    beta = 2 * pi * 0.2 * cfg.T;
    damp1 = sin(beta) / beta;
    damp2 = sin(2*beta) / (2*beta);
    Exact_Mean = 1 + 0.2 * sin(theta) * damp1;
    E_sin2 = 0.5 * (1 - cos(2 * theta) * damp2); % E[sin^2]
    E_sin_sq = (sin(theta) * damp1).^2;          % (E[sin])^2
    Exact_std = sqrt(0.04 * (E_sin2 - E_sin_sq));

    [num_Mean,num_Std] = gpc_SG4_1(cfg);
    L2_Error = L2Error(Exact_Mean, num_Mean);
    L2_errors_M(idx_M) = L2_Error;
    Linf_Error = LinfError(Exact_Mean, num_Mean);
    Linf_errors_M(idx_M) = Linf_Error;

    L2_Error_std = L2Error(Exact_std, num_Std);
    L2_errors_std(idx_M) = L2_Error_std;
    Linf_Error_std = LinfError(Exact_std, num_Std);
    Linf_errors_std(idx_M) = Linf_Error_std;

    fprintf('Grid Nx = %d, gPC Order M = %d\n', cfg.Nx, M);
    fprintf('L2 Error of Mean: %e\n', L2_Error);
    fprintf('Linf Error of Mean: %e\n', Linf_Error);
    fprintf('L2 Error of Std: %e\n', L2_Error_std);
    fprintf('Linf Error of Std: %e\n', Linf_Error_std);

    % 最后一个 M
    if idx_M == length(M_values)
        figure;
        plot(x_grid, Exact_Mean, 'b-', 'LineWidth', 1.5); hold on;
        plot(x_grid, num_Mean, 'r--', 'LineWidth', 1.5);
        legend('数值解期望', '期望解');
        title(['(t=' num2str(cfg.T) ')时的图像']);
        xlabel('x'); grid on;
    end
end

figure;
hold on;
semilogy(M_values, L2_errors_M, 'o-', 'LineWidth', 1.5,'MarkerSize', 8, 'DisplayName', 'Mean Error (L^2)');
semilogy(M_values, L2_errors_std, 's--','LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','Std Error (L^2)');
title('不同配置点个数下的 L^2 误差');
xlabel('gPC 阶 M');
ylabel('log(L^2 误差)');
set(gca, 'YScale', 'log');
ylim([1e-11, 1e-1]);                    
yticks([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]);
yticklabels({'10^{-12}', '10^{-11}', '10^{-10}', '10^{-9}', '10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'});
grid on;
legend('Location', 'Best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
hold off;

figure;
hold on;

semilogy(M_values, Linf_errors_M, 'o--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Mean Error (L^\infty)');
semilogy(M_values, Linf_errors_std, 's--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Variance Error (L^\infty)');
title('不同配置点个数下的 L_\infty 误差');
xlabel('gPC 阶 M');
ylabel('log(L_\infty 误差)');
set(gca, 'YScale', 'log');
ylim([1e-11, 1e-1]);                     
yticks([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]);
yticklabels({'10^{-12}', '10^{-11}', '10^{-10}', '10^{-9}', '10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'});
legend('Location', 'Best', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
hold off;




