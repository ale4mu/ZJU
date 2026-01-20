% 矩阵运算的辅助函数,批量计算jacobian矩阵和特征值
function [A_batch, alpha_vec] = Euler_Eigen_Batch(U_batch, gamma)
    % U_batch: (3, N_quad)
    % A_batch: (3, 3, N_quad)
    % alpha_vec: (1, N_quad)
    
    [~, Nq] = size(U_batch);
    
    rho = U_batch(1, :);
    m   = U_batch(2, :);
    E   = U_batch(3, :);
    eps = 1e-10;
    rho(rho < eps) = eps;

    u = m ./ rho;
    u2 = u.^2;
    p = (gamma - 1) * (E - 0.5 * rho .* u2);
    p(p < 1e-10) = 1e-10;
    
    c2 = gamma * p ./ rho;
    c2(c2 < 0) = 0; % 保护
    c = sqrt(c2);
    alpha_vec = abs(u) + c;
    
    % 批量构造 Jacobian 矩阵
    % A = [0, 1, 0; 
    %      0.5(g-3)u^2, (3-g)u, g-1;
    %      u(0.5(g-1)u^2 - H), H - (g-1)u^2, g*u]
    
    H = (E + p) ./ rho;
    gm1 = gamma - 1;
    A_batch = zeros(3, 3, Nq);
    
    % Row 1
    A_batch(1, 2, :) = 1;
    
    % Row 2
    A_batch(2, 1, :) = 0.5 * (gamma - 3) * u2;
    A_batch(2, 2, :) = (3 - gamma) * u;
    A_batch(2, 3, :) = gm1;
    
    % Row 3
    A_batch(3, 1, :) = u .* (0.5 * gm1 * u2 - H);
    A_batch(3, 2, :) = H - gm1 * u2;
    A_batch(3, 3, :) = gamma * u;
end