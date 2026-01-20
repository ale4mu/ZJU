% 欧拉方程的特征值
function [A, L, alpha, lambda] = Euler_Eigen(U, gamma)
    rho = U(1); m = U(2); E = U(3);
    u = m / rho;
    p = (gamma - 1) * (E - 0.5 * rho * u^2);
    c = sqrt(gamma * p / rho);
    alpha = abs(u) + c;
    
    % Jacobian A：αF/αU
    A = zeros(3,3);
    H = (E+p)/rho;
    gm1 = gamma - 1;
    u2 = u^2;
    
    A(1,2) = 1;
    A(2,1) = 0.5*(gamma-3)*u2; A(2,2) = (3-gamma)*u; A(2,3) = gm1;
    A(3,1) = u*(0.5*gm1*u2 - H); A(3,2) = H - gm1*u2; A(3,3) = gamma*u;

    % 左特征向量
    % 对应特征值: u-c, u, u+c
    L = zeros(3,3);
    
    b1 = 0.5 * gm1 * u2;
    b2 = gm1 * u;
    
    % u - c
    L(1,1) = (b1 + c*u) / (2*c^2);
    L(1,2) = -(b2 + c)  / (2*c^2);
    L(1,3) = gm1        / (2*c^2);
    
    % u 
    L(2,1) = 1 - b1 / c^2; 
    L(2,2) = b2 / c^2;
    L(2,3) = -gm1 / c^2;
    
    % u + c
    L(3,1) = (b1 - c*u) / (2*c^2);
    L(3,2) = -(b2 - c)  / (2*c^2);
    L(3,3) = gm1        / (2*c^2);
    
    % 特征值
    lambda = diag([u-c, u, u+c]);
    
    % check = norm(L*A - diag(lambda)*L);
    % if check > 1e-10
    %     error('L not consistent with A');
    % end
end