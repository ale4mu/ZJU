function [A_batch, alpha_vec] = Euler_Eigen_Batch_2D(U_batch, gamma)
    % U_batch: (4, N_quad) -> [rho, rho*u, rho*v, E]
    [nVar, Nq] = size(U_batch);
    
    rho = U_batch(1, :);
    m_x = U_batch(2, :); % rho * u
    m_y = U_batch(3, :); % rho * v
    E   = U_batch(4, :);
    
    epsilon = 1e-2;
    mask_rho = rho < epsilon;
    rho(mask_rho) = epsilon;
    % 如果密度修正了，动量也要清零防止速度过大
    m_x(mask_rho) = 0;
    m_y(mask_rho) = 0;
    
    u = m_x ./ rho;
    v = m_y ./ rho; 
    q2 = u.^2 + v.^2; 
    p = (gamma - 1) * (E - 0.5 * rho .* q2);
    
    mask_p = p < epsilon;
    p(mask_p) = epsilon;
    % 如果压力修正了，能量也要修正
    E(mask_p) = p(mask_p)/(gamma-1) + 0.5 * rho(mask_p) .* q2(mask_p);
    
    c2 = gamma * p ./ rho;
    c = sqrt(c2);
    alpha_vec = abs(u) + c; % 特征值只与 u 有关 (x方向)
    
    % 4x4 Jacobian 矩阵
    H = (E + p) ./ rho;
    gm1 = gamma - 1;
    
    A_batch = zeros(4, 4, Nq);
    
    A_batch(1, 2, :) = 1;
    A_batch(2, 1, :) = 0.5 * gm1 * q2 - u.^2;
    A_batch(2, 2, :) = (3 - gamma) * u;
    A_batch(2, 3, :) = -gm1 * v;
    A_batch(2, 4, :) = gm1;
    A_batch(3, 1, :) = -u .* v;
    A_batch(3, 2, :) = v;
    A_batch(3, 3, :) = u;
    % A_batch(3, 4, :) = 0;

    A_batch(4, 1, :) = u .* (0.5 * gm1 * q2 - H);
    A_batch(4, 2, :) = H - gm1 * u.^2;
    A_batch(4, 3, :) = -gm1 * u .* v;
    A_batch(4, 4, :) = gamma * u;
end