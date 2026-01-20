function [Mean_Rho,Mean_Rho2] = colloaction41(cfg,initial_func)
    % cfg: PDE的参数
    cfg.Nx = 320;
    % 获取 Gauss-Legendre 积分点和权重
    num_points = 20;
    [nodes, weights] = lgwt(num_points, -1, 1);
    Results = zeros(num_points, cfg.Nx);
    
    % 初始化统计量
    Mean_Rho = zeros(1, cfg.Nx);
    Mean_Rho2 = zeros(1, cfg.Nx); % E[rho^2]
    %parpool('local',2);
    parfor j = 1:num_points
        xi = nodes(j);
        % 用当前配置点xi求解确定性方程
        current_init_func = @(x) initial_func(x, cfg, xi); 
        Results(j, :) = solve_SC(cfg, current_init_func);
    end

    for j = 1:num_points
        w_j = weights(j) / 2.0; % 因为 xi ~ U[-1, 1]
        res = Results(j, :);
        Mean_Rho = Mean_Rho + res * w_j;
        Mean_Rho2 = Mean_Rho2 + (res.^2) * w_j;
    end
end
    