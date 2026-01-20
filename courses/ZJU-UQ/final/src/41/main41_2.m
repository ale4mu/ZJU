% 用于验证WENO的5阶收敛阶
clear;clc;close all;

cfg.gamma = 1.4;       % 气体常数
cfg.L = 1.0;           % 区域长度
cfg.T = 0.2;
cfg.CFL = 0.6;
NX = [10,20,40,80,160];
M = 10;
cfg.M = M;
L2_errors_M = zeros(size(NX));
Linf_errors_M = zeros(size(NX));
L2_errors_VAR = zeros(size(NX));
Linf_errors_VAR = zeros(size(NX));

initial_func = @initial_func4_1;

for idx_NX = 1:length(NX)
    cfg.Nx = NX(idx_NX);
    
    [num_Mean_rho,num_Std_rho] = gpc_SG4_1(cfg);
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

    L2_Error = L2Error(num_Mean_rho, Exact_Mean);
    L2_errors_M(idx_NX) = L2_Error;
    Linf_Error = LinfError(num_Mean_rho, Exact_Mean);
    Linf_errors_M(idx_NX) = Linf_Error;

    L2_Error_Var = L2Error(num_Std_rho, Exact_std);
    L2_errors_VAR(idx_NX) = L2_Error_Var;
    Linf_Error_Var = LinfError(num_Std_rho, Exact_std);
    Linf_errors_VAR(idx_NX) = Linf_Error_Var;
end

fprintf('%-10s | %-18.2s | %-12.4s | %-10.2s | %-18.4s | %-12.2s | %-10.4s | %-18.2s | %-12.4s\n', ...
    'NX', 'L2(Mean)', 'L2 Order', ...
    'Linf(Mean)', 'Linf Order', 'L2(Std)', ...
    'L2 Order', 'Linf(Std)', 'Linf Order');

fprintf('%-10d | %-18.2e | %-12.4f | %-10.2e | %-18.4f | %-12.2e | %-10.4f | %-18.2e | %-12.4f\n', ...
    NX(1), L2_errors_M(1), 0, Linf_errors_M(1), 0, ...
    L2_errors_VAR(1), 0, Linf_errors_VAR(1), 0);

for i = 2:length(NX)
    l2_order_mean   = order(L2_errors_M(i-1),   L2_errors_M(i),   NX(i-1), NX(i));
    linf_order_mean = order(Linf_errors_M(i-1), Linf_errors_M(i), NX(i-1), NX(i));
    l2_order_std    = order(L2_errors_VAR(i-1), L2_errors_VAR(i), NX(i-1), NX(i));
    linf_order_std  = order(Linf_errors_VAR(i-1), Linf_errors_VAR(i), NX(i-1), NX(i));

    fprintf('\n%s\n', repmat('-', 1, 135));

    fprintf('%-10d | %-18.2e | %-12.4f | %-10.2e | %-18.4f | %-12.2e | %-10.4f | %-18.2e | %-12.4f\n', ...
        NX(i), ...
        L2_errors_M(i),     l2_order_mean, ...
        Linf_errors_M(i),   linf_order_mean, ...
        L2_errors_VAR(i),   l2_order_std, ...
        Linf_errors_VAR(i), linf_order_std);
end
