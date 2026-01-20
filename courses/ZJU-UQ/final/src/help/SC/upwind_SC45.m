function RHS = upwind_SC45(U, dx, gamma)
    % 1阶迎风格式
    % U: (4, Nx)
    [nVar, Nx] = size(U);
    rho = U(1, :); m_x = U(2, :); m_y = U(3, :); E = U(4, :);
    eps = 1e-6;
    rho(rho < eps) = eps;
    
    u = m_x ./ rho;
    v = m_y ./ rho;
    q2 = u.^2 + v.^2;
    p = (gamma-1)*(E - 0.5*rho.*q2);
    p(p < eps) = eps;
    
    c = sqrt(gamma*p./rho);
    alpha = abs(u) + c; 

    F = zeros(4, Nx);
    F(1, :) = m_x;
    F(2, :) = m_x.*u + p;
    F(3, :) = m_x.*v; % rho*u*v
    F(4, :) = u.*(E+p);
    
    % F_{i+1/2}
    F_hat = zeros(4, Nx+1);
    
    % 内部界面 i=2:Nx (对应 U(:, i-1) 和 U(:, i))
    for i = 2:Nx
        U_L = U(:, i-1);
        U_R = U(:, i);
        F_L = F(:, i-1);
        F_R = F(:, i);
        
        a_max = max(alpha(i-1), alpha(i));
        F_hat(:, i) = 0.5 * (F_L + F_R - a_max * (U_R - U_L));
    end
    
    % 边界处理
    % F_{1/2} (左边界) -> 使用 U(:, 1)
    % F_{N+1/2} (右边界) -> 使用 U(:, Nx)
    % 左边界通量就是物理通量 F(:, 1)
    F_hat(:, 1) = F(:, 1); 
    F_hat(:, Nx+1) = F(:, Nx);
    
    RHS = - (F_hat(:, 2:end) - F_hat(:, 1:end-1)) / dx;
end