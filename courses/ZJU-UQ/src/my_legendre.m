function [nodes, weights] = my_legendre(M)
    % GET_LEGENDRE_GAUSS 生成 Gauss-Legendre 积分节点和权重
    % 积分区间 [-1, 1], 权函数 w(x) = 1
    
    if M == 1
        nodes = 0;
        weights = 2; % 积分区间长度
        return;
    end

    % 1. 构建 Jacobi 矩阵
    % Legendre 多项式递推关系: 
    % (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
    % 归一化后的 Jacobi 矩阵次对角线系数 beta_n = n / sqrt(4n^2 - 1)
    
    n = 1:(M-1);
    beta = n ./ sqrt(4 * n.^2 - 1);
    
    % 对称三对角矩阵，对角线为0
    J = diag(beta, 1) + diag(beta, -1);
    
    % 2. 计算特征值和特征向量
    [V, D] = eig(J);
    
    % 节点即为特征值
    nodes = diag(D);
    
    % 3. 计算权重
    % Golub-Welsch 公式: w_j = 2 * (V(1,j))^2
    % 系数 2 来自于 int_{-1}^{1} 1 dx = 2 (Legendre权函数的0阶矩)
    weights = 2 * (V(1, :)'.^2);
    
    % 4. 排序输出
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
    
    % 确保为列向量
    nodes = nodes(:);
    weights = weights(:);
end