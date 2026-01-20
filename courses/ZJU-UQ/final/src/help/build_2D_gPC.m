function [Phi_val, nodes_2d, weights_2d, norm_sq] = build_2D_gPC(M, num_points_1d)
    % 构建2维gPC多项式
    % M: 总阶数
    % num_points_1d: 每个维度的积分点数
    
    % 生成 1D 积分点和权重
    [nodes_1d, weights_1d] = lgwt(num_points_1d, -1, 1);
    
    % 生成 2D 张量积积分点
    [XI1, XI2] = meshgrid(nodes_1d, nodes_1d);
    nodes_2d = [XI1(:), XI2(:)]'; % (2, N_quad_2d)
    
    [W1, W2] = meshgrid(weights_1d, weights_1d);
    weights_2d = (W1 .* W2); 
    weights_2d = weights_2d(:); % (N_quad_2d, 1)
    
    N_quad_2d = length(weights_2d);
    
    % 构建多重索引
    % 找出所有满足 i + j <= M 的 (i, j)
    idx_set = [];
    for i = 0:M
        for j = 0:M
            if i + j <= M
                idx_set = [idx_set; i, j];
            end
        end
    end
    nModes = size(idx_set, 1); % 基函数总数
    
    % 计算基函数值 Phi_val (nModes, N_quad_2d)
    % 和归一化因子 norm_sq (nModes, 1)
    Phi_val = zeros(nModes, N_quad_2d);
    norm_sq = zeros(nModes, 1);
    
    % 预计算 1D Legendre 值
    P_val_1 = zeros(M+1, N_quad_2d);
    P_val_2 = zeros(M+1, N_quad_2d);
    for deg = 0:M
        P_val_1(deg+1, :) = my_legendre(deg, nodes_2d(1, :));
        P_val_2(deg+1, :) = my_legendre(deg, nodes_2d(2, :));
    end
    
    for k = 1:nModes
        deg1 = idx_set(k, 1);
        deg2 = idx_set(k, 2);
        
        % Phi_k(xi) = P_deg1(xi1) * P_deg2(xi2)
        Phi_val(k, :) = P_val_1(deg1+1, :) .* P_val_2(deg2+1, :);
        
        % Norm_k = Norm_deg1 * Norm_deg2
        % Norm_1d = 2 / (2*n + 1)
        n1 = 2 / (2*deg1 + 1);
        n2 = 2 / (2*deg2 + 1);
        norm_sq(k) = n1 * n2;
    end
    fprintf('2D gPC Built: Total Degree M=%d, nModes=%d, N_quad=%d\n', M, nModes, N_quad_2d);
end