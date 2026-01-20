function [nodes, weights] = hermite(M)
   
    % 针对权函数 w(x) = 1/sqrt(2*pi) * exp(-x^2/2)
    % 积分区间 (-inf, inf)
    
    if M == 1
        nodes = 0;
        weights = 1;
        return;
    end
    
    %% 直接构建概率型 Jacobi 矩阵
    % 概率型递推系数 sqrt(n)，n = 1, 2, ..., M-1
    % 对比：物理型是 sqrt(n/2)
    beta = sqrt(1:M-1); 
    
    % 构建对称三对角矩阵
    J = diag(beta, 1) + diag(beta, -1);
    
    % 计算特征值和特征向量
    [V, D] = eig(J);
    
    % 节点即为特征值
    nodes = diag(D);
    
    % 权重计算
    % 标准 Golub-Welsch 公式: w_j = mu_0 * (v_{j,1})^2
    % 对于标准正态分布，0阶矩 mu_0 = integral(PDF) = 1
    % 所以不需要像物理型那样乘以 sqrt(pi)
    weights = (V(1, :).^2)';
    
    %% 整理输出
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
    
    % 确保为列向量
    nodes = nodes(:);
    weights = weights(:);
end