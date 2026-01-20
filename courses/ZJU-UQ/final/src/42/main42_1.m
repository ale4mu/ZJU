% 用于验证WENO的5阶收敛阶
clear;clc;close all;

cfg.gamma = 1.4;       % 气体常数
cfg.L = 1.0;           % 区域长度
cfg.T = 0.2;
cfg.CFL = 0.6;
NX = [10,20,40,80];
cfg.M = 6;
cfg_ref = cfg;
cfg_ref.Nx = 320;

[Mean_Ref, Std_Ref] = collocation42(cfg_ref);


L2_errors_M = zeros(size(NX));
Linf_errors_M = zeros(size(NX));
L2_errors_VAR = zeros(size(NX));
Linf_errors_VAR = zeros(size(NX));

for idx_NX = 1:length(NX)
    nx = NX(idx_NX);
    cfg.Nx = nx;

    [num_Mean_rho,num_Std_rho] = gPC_SG4_2(cfg);
    dx = cfg.L / cfg.Nx;
    x_grid = (0.5:cfg.Nx-0.5) * dx;
    step = cfg_ref.Nx / nx;

    Mean_Ref_Down = zeros(1, nx);
    Std_Ref_Down  = zeros(1, nx);
    for i = 1:nx
        % 粗网格单元 i 对应细网格的范围
        idx_start = (i-1)*step + 1;
        idx_end   = i*step;
        % 取平均 (有限体积守恒律)
        Mean_Ref_Down(i) = mean(Mean_Ref(idx_start:idx_end));
        Std_Ref_Down(i)  = sqrt(mean(Std_Ref(idx_start:idx_end).^2));
        %Std_Ref_Down(i)  = mean(Std_Ref(idx_start:idx_end)); 
    end

    L2_Error = L2Error(num_Mean_rho, Mean_Ref_Down);
    L2_errors_M(idx_NX) = L2_Error;
    Linf_Error = LinfError(num_Mean_rho, Mean_Ref_Down);
    Linf_errors_M(idx_NX) = Linf_Error;

    L2_Error_Var = L2Error(num_Std_rho, Std_Ref_Down);
    L2_errors_VAR(idx_NX) = L2_Error_Var;
    Linf_Error_Var = LinfError(num_Std_rho, Std_Ref_Down);
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
