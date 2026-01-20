% 用于验证谱收敛
clear;clc;close all;

cfg.gamma = 1.4;       % 气体常数（绝热率）
cfg.L = 1.0;           % 区域长度
cfg.T = 0.2;
cfg.CFL = 0.6;
cfg.Nx = 80;

M_values = [0,1, 2, 3, 4]; 
L2_errors_M = zeros(size(M_values));
Linf_errors_M = zeros(size(M_values));
L2_errors_std = zeros(size(M_values));
Linf_errors_std = zeros(size(M_values));


for idx_M = 1:length(M_values)
    M = M_values(idx_M);
    cfg.M = M;
    
    dx = cfg.L / cfg.Nx;
    x_grid = (0.5:cfg.Nx-0.5) * dx;

    [num_Mean,num_Std] = gPC_SG4_2(cfg);
    [Mean_Ref, Std_Ref] = collocation42(cfg);
    L2_Error = L2Error(Mean_Ref, num_Mean);
    L2_errors_M(idx_M) = L2_Error;
    Linf_Error = LinfError(Mean_Ref, num_Mean);
    Linf_errors_M(idx_M) = Linf_Error;

    L2_Error_std = L2Error(Std_Ref, num_Std);
    L2_errors_std(idx_M) = L2_Error_std;
    Linf_Error_std = LinfError(Std_Ref, num_Std);
    Linf_errors_std(idx_M) = Linf_Error_std;

    fprintf('Grid Nx = %d, gPC Order M = %d\n', cfg.Nx, M);
    fprintf('L2 Error of Mean: %e\n', L2_Error);
    fprintf('Linf Error of Mean: %e\n', Linf_Error);
    fprintf('L2 Error of Std: %e\n', L2_Error_std);
    fprintf('Linf Error of Std: %e\n', Linf_Error_std);

end

figure;
hold on;
semilogy(M_values, L2_errors_M, 'o-', 'LineWidth', 1.5,'MarkerSize', 8, 'DisplayName', 'Mean Error (L2)');
semilogy(M_values, L2_errors_std, 's--','LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','Std Error (L2)');
xlabel('gPC 阶 M');
ylabel('log(L2误差)');
set(gca, 'YScale', 'log');
ylim([1e-8, 1e-1]);                    
yticks([1e-8, 1e-7, 1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]);
yticklabels({'10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'});
grid on;
legend('Location', 'Best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
hold off;

figure;
hold on;

semilogy(M_values, Linf_errors_M, 'o--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Mean Error (L_\infty)');
semilogy(M_values, Linf_errors_std, 's--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Variance Error (L_\infty)');
xlabel('gPC 阶 M');
ylabel('log(L_\infty 误差)');
set(gca, 'YScale', 'log');
ylim([1e-8, 1e-1]);                     
yticks([1e-8, 1e-7, 1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]);
yticklabels({'10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}'});
legend('Location', 'Best', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12);
hold off;




