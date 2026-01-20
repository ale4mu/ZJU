function [U_L, U_R] = WENO_gPC_2D(U_hat, Phi_val, gamma, nodes, weights, norm_sq)
    % Example 4.5
    % 边界条件：零梯度外推?
    
    [nVar, nModes, Nx] = size(U_hat);
    n_L = 3; 
    n_R = 4;
    
    U_ext = zeros(nVar, nModes, Nx + n_L + n_R);
    
    % 中间
    U_ext(:, :, n_L+1 : n_L+Nx) = U_hat;
    
    % 左边界 (复制第一个单元)
    for g = 1:n_L
        U_ext(:, :, g) = U_hat(:, :, 1);
    end
    
    % 右边界 (复制最后一个单元)
    for g = 1:n_R
        U_ext(:, :, n_L + Nx + g) = U_hat(:, :, Nx);
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
                % 正通量重构 U_L (界面 i+1/2 左侧)
                stencil_pos = squeeze(U_ext(v, k, idx_curr-2 : idx_curr+2));
                u_l(v, k) = WENO_pos(stencil_pos);
                
                % 负通量重构 U_R (界面 i+1/2 右侧)
                stencil_neg = squeeze(U_ext(v, k, idx_curr : idx_curr+4));
                u_r(v, k) = WENO_neg(stencil_neg);
            end
        end
        
        U_L(:, :, i+1) = Limiter(u_l, U_avg_curr, Phi_val, gamma);
        U_R(:, :, i+1) = Limiter(u_r, U_avg_next, Phi_val, gamma);
    end
    
    % 左边界是常数外推，所以 U_L(1) = U_hat(1)
    U_L(:, :, 1) = U_hat(:, :, 1);
    
    % U_R(1) 是第一个单元的左重构值
    idx_1 = 1 + n_L;
    U_avg_1 = U_ext(:, :, idx_1);
    u_r_1 = zeros(nVar, nModes);
    for k = 1:nModes
        for v = 1:nVar
            stencil_neg_1 = squeeze(U_ext(v, k, idx_1-1 : idx_1+3));
            u_r_1(v, k) = WENO_neg(stencil_neg_1);
        end
    end
    U_R(:, :, 1) = Limiter(u_r_1, U_avg_1, Phi_val, gamma);

end