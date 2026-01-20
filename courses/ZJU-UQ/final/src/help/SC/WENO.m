%% L(U)
% 输入：状态 U
% 输出：dU/dt 
function RHS = WENO(U, dx, gamma)
    [nVar, Nx] = size(U);
    
    % U = (rho,rho * u, E)
    rho = U(1,:);
    m   = U(2,:); % rho * u
    E   = U(3,:);
    u   = m ./ rho; % 逐元素除法
    % E = rho * e + 0.5 * rho * u^2
    % p = (gamma - 1) * rho * e
    p   = (gamma - 1) * (E - 0.5 * rho .* u.^2); 
    
    % F(U) = (rho * u, rho * u^2 + p, uE + up)
    F = zeros(size(U));
    F(1,:) = m;
    F(2,:) = m.*u + p;
    F(3,:) = (E + p).*u;
    
    % Lax-Friedrichs 分裂
    % alpha：Jacobbi矩阵的最大特征值的绝对值（？）
    c = sqrt(gamma * p ./ rho);
    alpha = max(abs(u) + c); % 全局Lax-Friedrichs 分裂，如果效果不好可以使用局部Lax-Friedrichs 分裂
    
    F_pos = 0.5 * (F + alpha * U); % 正通量 F+
    F_neg = 0.5 * (F - alpha * U); % 负通量 F- 
    
    % 计算数值通量 F_hat_{j+1/2}
    F_hat = zeros(nVar, Nx+1);
    
    % 处理周期性边界条件
    % 扩展数组: [最后3个, 原始数据, 最前3个] -> 用于索引 i-3 到 i+3
    idx_2 = [Nx-2,Nx-1, Nx, 1:Nx, 1, 2,3]; 
    F_pos_ext = F_pos(:, idx_2);
    F_neg_ext = F_neg(:, idx_2);
    
    % 循环计算每一个界面 j+1/2 的通量
    % 界面 i 对应原始网格 i 和 i+1 之间
    for i = 1:Nx
        idx = i + 3; % 前面加了3个
        % 对每个U(rho, rho*u, E)分别做WENO
        for k = 1:nVar
            % F+: v_{j-2}, ..., v_{j+2} -> 得到界面左侧值
            v_pos = F_pos_ext(k, idx-2 : idx+2);
            f_plus = WENO_pos(v_pos); % f+
            
            % F-: v_{j+3}, ..., v_{j-1} (镜像对称) -> 得到界面右侧值
            v_neg = F_neg_ext(k, idx-1 : idx+3);
            f_minus = WENO_neg(v_neg); % f-

            
            F_hat(k, i+1) = f_plus + f_minus;
        end
    end

    % F_hat(:, 1) （即F_{1/2}） 应该等于 F_{N+1/2} (即 F_hat(:, Nx+1))
    F_hat(:, 1) = F_hat(:, Nx+1);
    
    % 计算残差 RHS = - (F_{j+1/2} - F_{j-1/2}) / dx
    % 对应网格 j，右界面是 j+1 (索引i+1)，左界面是 j (索引i)
    RHS = - (F_hat(:, 2:Nx+1) - F_hat(:, 1:Nx)) / dx;
end