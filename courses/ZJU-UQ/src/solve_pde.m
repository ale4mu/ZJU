function [err_mean, err_var] = solve_pde(N, M, T)
    % 1. 获取 Legendre 节点和权重 ([-1, 1] 区间)
    [z, w] = my_legendre(M);
    
    % C++: w_j = weights[j] / 2.0; 
    % 原因：z ~ U[-1, 1]，PDF = 1/2。期望积分 E[f] = int f(z)*0.5 dz
    w = w / 2.0;
    
    % 空间网格设置 (假设区间 [0, 2*pi])
    L = 2 * pi;
    dx = L / N;
    x = (0:N-1)' * dx; % 对应C++ vector索引 0 到 N-1
    
    % 初始化统计量
    y_mean_num = zeros(N, 1);
    y_sq_mean_num = zeros(N, 1);
    
    % 2. 随机配点循环
    for j = 1:M
        % C++: pde(1 + 0.1 * nodes[j], ...)
        % 随机参数 a = 1 + 0.1 * z
        a = 1.0 + 0.1 * z(j);
        
        % 调用数值格式求解器 scheme3
        % 注意：这里假设 scheme3 返回的是 T 时刻的空间分布向量 (Nx1)
        % 你需要确保你的 scheme3 接受 (a, N, dx, T) 或者类似的参数
        [~,u_full] = scheme3(N,T,a,1); 

        u_val = u_full(1:N);
        u_val = u_val(:); 
        
        % 累加期望和二阶矩
        y_mean_num = y_mean_num + u_val * w(j);
        y_sq_mean_num = y_sq_mean_num + (u_val.^2) * w(j);
    end
    
    % 计算数值解方差
    y_var_num = y_sq_mean_num - y_mean_num.^2;
    
    % 3. 计算解析解 (Exact Solution)
    % C++ Lambda: pde_exact
    % exact_mean = sin(x + T) * sin(0.1 * T) / (0.1 * T);
    term_mean = sin(0.1 * T) / (0.1 * T);
    exact_mean = sin(x + T) * term_mean;
    
    % exact_var = 0.5 * (1 - cos(2 * (x + T)) * sin(0.2 * T) / (0.2 * T)) - exact_mean^2;
    term_var = sin(0.2 * T) / (0.2 * T);
    exact_var = 0.5 * (1 - cos(2 * (x + T)) * term_var) - exact_mean.^2;
    
    % 4. 计算 L2 误差
    % err = sqrt( sum(diff^2) * dx )
    diff_mean = y_mean_num - exact_mean;
    diff_var = y_var_num - exact_var;
    
    err_mean = sqrt(sum(diff_mean.^2) * dx);
    err_var = sqrt(sum(diff_var.^2) * dx);
end