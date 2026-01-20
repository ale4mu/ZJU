function [L, Lambda] = Euler_Eigen_2D(U, gamma)
    % 计算 4 变量 Euler 方程的左特征向量矩阵 L 和特征值对角阵 Lambda
    % U = [rho, m_x, m_y, E]
    
    rho = U(1); m = U(2); my = U(3); E = U(4);
    
    if rho < 1e-4, rho = 1e-4; end
    
    u = m / rho;
    v = my / rho;
    q2 = u^2 + v^2;
    p = (gamma - 1) * (E - 0.5 * rho * q2);
    
    if p < 1e-4, p = 1e-4; end
    
    c = sqrt(gamma * p / rho);
    
    % 特征值: u-c, u, u, u+c
    Lambda = diag([u-c, u, u, u+c]);
    
    % 左特征向量矩阵 L (行向量是左特征向量)
    % 对应变量顺序: rho, rho*u, rho*v, E
    
    gm1 = gamma - 1;
    L = zeros(4, 4);
    
    % u - c 
    b1 = 0.5 * gm1 * q2;
    b2 = gm1 * u;
    b3 = gm1 * v;
    
    L(1, 1) = (b1 + c*u) / (2*c^2);
    L(1, 2) = -(b2 + c)  / (2*c^2);
    L(1, 3) = -b3        / (2*c^2);
    L(1, 4) = gm1        / (2*c^2);
    
    % u (熵波 - 密度变化)
    % 对应 eigenvector 使得 dP = 0
    L(2, 1) = 1 - b1/c^2;
    L(2, 2) = b2/c^2;
    L(2, 3) = b3/c^2;
    L(2, 4) = -gm1/c^2;
    
    % u (剪切波 - v 变化)
    % 这是一个线性简并场，只涉及 v 的变化，不影响 p 和 u
    % 构造一个正交于其他向量的行
    L(3, 1) = -v;
    L(3, 2) = 0;
    L(3, 3) = 1;
    L(3, 4) = 0;
    
    % u + c (声波)
    L(4, 1) = (b1 - c*u) / (2*c^2);
    L(4, 2) = -(b2 - c)  / (2*c^2);
    L(4, 3) = -b3        / (2*c^2);
    L(4, 4) = gm1        / (2*c^2);
end