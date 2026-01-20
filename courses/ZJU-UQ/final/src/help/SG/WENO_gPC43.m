function [U_L, U_R] = WENO_gPC43(U_hat, Phi_val, gamma, current_time, nodes, weights, norm_sq)
    % U_hat: (nVar, nModes, Nx)
    % Phi_val: (nModes, N_quad)
    % nodes: (N_quad, 1) 积分点
    % weights: (N_quad, 1) 权重
    % norm_sq: (nModes, 1) 归一化因子
    [nVar, nModes, Nx] = size(U_hat);
    
    % u(0,t,xi) = 0.02*sin(2*pi*w(xi)*t)
    % w(xi) = 1 + 0.1*xi
    w_xi = 1 + 0.1 * nodes; 
    u_phys_boundary = 0.02 * sin(2 * pi * w_xi * current_time); 
    
    % 计算物理空间的能量 E_phys
    % E = p/(gamma-1) + 0.5 * rho * u^2
    rho_phys = 1.0;
    p_phys   = 0.6;
    E_phys_boundary = p_phys / (gamma - 1) + 0.5 * rho_phys * (u_phys_boundary.^2);
    
    % 投影得到 gPC 系数
    % U_BC_coeffs: (3, nModes) -> [rho_hat; m_hat; E_hat]
    U_BC_coeffs = zeros(nVar, nModes);
    
    % rho 的系数 (常数 1)
    U_BC_coeffs(1, 1) = 1.0; % (其他高阶模态为0，初始化已默认为0)
    
    % m = rho*u 的系数 (因为 rho=1, 所以 m=u)
    % E 的系数
    for k = 1:nModes
        basis = Phi_val(k, :); 
        
        vec_basis = basis(:);
        vec_u     = u_phys_boundary(:);
        vec_E     = E_phys_boundary(:);
        vec_w     = weights(:);
        
        % 投影公式: <f, phi_k> / <phi_k, phi_k>
        % sum( f(q) * w(q) * phi_k(q) )
        temp_m = sum(vec_u .* vec_w .* vec_basis);
        U_BC_coeffs(2, k) = temp_m / norm_sq(k);
        
        % 能量 E
        temp_E = sum(vec_E .* vec_w .* vec_basis);
        U_BC_coeffs(3, k) = temp_E / norm_sq(k);
    end

    % 扩展数组 
    n_L = 3;
    n_R = 4;
    U_ext = zeros(nVar, nModes, Nx + n_L + n_R);
    U_ext(:, :, n_L+1 : n_L+Nx) = U_hat;
    
    % 左边界
    for g = 1:n_L
        U_ext(:, :, g) = U_BC_coeffs;
    end
    
    % 右边界，简单复制
    U_last = U_hat(:, :, Nx);
    for g = 1:n_R
        U_ext(:, :, n_L + Nx + g) = U_last;
    end

    % WENO 重构
    U_L = zeros(nVar, nModes, Nx+1);
    U_R = zeros(nVar, nModes, Nx+1);
    
    for i = 1:Nx
        idx_curr = i + n_L; 
        idx_next = idx_curr + 1;
        U_avg_curr = U_ext(:, :, idx_curr); 
        U_avg_next = U_ext(:, :, idx_next);
        u_l = zeros(nVar, nModes);
        u_r = zeros(nVar, nModes);
        for k = 1:nModes
            for v = 1:nVar
                stencil_pos = squeeze(U_ext(v, k, idx_curr-2 : idx_curr+2));
                u_l(v, k) = WENO_pos(stencil_pos); 
                stencil_neg = squeeze(U_ext(v, k, idx_curr : idx_curr+4));
                u_r(v, k) = WENO_neg(stencil_neg);
            end
        end
        % Examle 4.3需要使用限制器，不然数值会崩溃
        U_L(:, :, i+1) = Limiter(u_l, U_avg_curr, Phi_val, gamma);
        U_R(:, :, i+1) = Limiter(u_r, U_avg_next, Phi_val, gamma);
    end

    % 最左边界
    U_L(:, :, 1) = U_BC_coeffs;
    idx_1 = 1 + n_L;
    U_avg_1 = U_ext(:, :, idx_1);

    u_r_1 = zeros(nVar, nModes);
    for k = 1:nModes
        for v = 1:nVar
            stencil_neg_1 = squeeze(U_ext(v, k, idx_1-1 : idx_1+3));
            u_r_1(v, k) = WENO_neg(stencil_neg_1);
        end
    end
    % 使用限制器
    U_R(:, :, 1) = Limiter(u_r_1, U_avg_1, Phi_val, gamma);
end