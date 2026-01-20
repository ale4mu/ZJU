% Example 4.3
function RHS = upwind_SC43(U, t, w_freq, dx, gamma)
    [~, Nx] = size(U);
    
    % 左边界: u(0,t) = 0.02 * sin(2*pi*w*t)
    u_bc = 0.02 * sin(2 * pi * w_freq * t);
    rho_bc = 1.0;
    p_bc = 0.6;
    E_bc = p_bc/(gamma-1) + 0.5*rho_bc*u_bc^2;
    U_left_BC = [rho_bc; rho_bc*u_bc; E_bc];
    
    n_L = 3;
    n_R = 4;
    U_ext = zeros(3, Nx + n_L + n_R);
    U_ext(:, n_L+1:n_L + Nx) = U;
  
    for g = 1:n_L
        U_ext(:, g) = U_left_BC; 
    end
    for g = 1:n_R
        U_ext(:, n_L + Nx + g) = U(:, end); 
    end
    
    F_hat = zeros(3, Nx+1);
    rho = U_ext(1,:); 
    m = U_ext(2,:); 
    E = U_ext(3,:);
    u = m./rho; 
    p = (gamma-1)*(E - 0.5*rho.*u.^2);
    c = sqrt(gamma*abs(p)./rho);
    alpha = max(abs(u)+c);
    
    F = zeros(size(U_ext));
    F(1,:) = m;
    F(2,:) = m.*u + p;
    F(3,:) = (E+p).*u;
    
    F_pos = 0.5*(F + alpha*U_ext);
    F_neg = 0.5*(F - alpha*U_ext);
    
    for i = 1:Nx+1
        idx = i + n_L; % 对应界面 i-1/2 在 U_ext 中的位置 
        for k = 1:3
            % 使用 1 阶迎风格式 
            % F+ : 取 idx-1
            fp = F_pos(k, idx-1);
            % F- : 取 idx
            fn = F_neg(k, idx);
            F_hat(k, i) = fp + fn;
        end
    end
    
    RHS = - (F_hat(:, 2:end) - F_hat(:, 1:end-1)) / dx;
end