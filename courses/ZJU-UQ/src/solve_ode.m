%% 4. 求解ODE的函数
function [mean_error, var_error] = solve_ode(N, M)
    % 用随机配点法 + Runge-Kutta法求解ODE
    % N: 时间离散度
    % M: 随机配置点个数
    
    % 时间步长
    T = 1;
    beta = 1;
    dt = T / N;
    
    % 获取概率型Hermite多项式的节点和权重
    [nodes, weights] = hermite(M);
    
    % 初始化统计量
    mean_y = zeros(1, N+1);
    mean_y2 = zeros(1, N+1);
    
    % 遍历所有配置点
    for j = 1:M
        alpha = nodes(j);
        w = weights(j);

        % 求解确定性ODE
        [~,y] = RK_three(N,0,T,beta,alpha);
        
        % 累积统计量
        mean_y = mean_y + w * y;
        mean_y2 = mean_y2 + w * (y.^2);
    end
    
    % 计算方差
    var_y = mean_y2 - mean_y.^2;
    
    % 时间点
    t_points = 0:dt:T;
    
    % 计算精确解
    [exact_mean, exact_var] = exact_solution(t_points);
    
    % 计算误差（使用L∞范数，即最大绝对误差）
    mean_error = max(abs(mean_y - exact_mean));
    var_error = max(abs(var_y - exact_var));
end