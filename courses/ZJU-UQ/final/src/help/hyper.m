% 检查是否保双曲性，即特征值是否都是实数
function hyper(U_hat_slice, gamma, Phi_val, weights, norm_sq)
    % U_hat_slice: 某个空间点上的系数 (3, M+1)
    % 验证该点处的特征值是否为实数
    [nVar, nModes] = size(U_hat_slice);
    N_quad = length(weights);  

    % 验证是否处在容许状态集中
    [is_valid, ~, ~] = check_Euler(U_hat_slice, Phi_val, gamma);
    if ~is_valid
        error('Invalid solution!');
    end
    U_phys = U_hat_slice * Phi_val; % (3, N_quad)
    
    % 计算 A0 和 A1
    % A0_ij = 积分 A0(U) phi_k phi_l du(xi)
    A0_hat = zeros(nVar*nModes, nVar*nModes);
    A1_hat = zeros(nVar*nModes, nVar*nModes);
    
    for q = 1:N_quad
        u_val = U_phys(:, q);
        w = weights(q);
        
        % 计算物理空间的 A0, A1
        % A0 = L^T * L, A1 = L^T * Lambda * L
        [~, L, ~,Lambda] = Euler_Eigen(u_val, gamma); 
        
        A0_phys = L' * L;
        A1_phys = L' * Lambda * L;
        
        % 累加到 Galerkin 矩阵
        for k = 1:nModes
            for l = 1:nModes
                % 块索引
                r_idx = (k-1)*nVar + (1:nVar);
                c_idx = (l-1)*nVar + (1:nVar);
                
                basis_prod = Phi_val(k, q) * Phi_val(l, q);
                
                A0_hat(r_idx, c_idx) = A0_hat(r_idx, c_idx) + w * A0_phys * basis_prod;
                A1_hat(r_idx, c_idx) = A1_hat(r_idx, c_idx) + w * A1_phys * basis_prod;
            end
        end
    end
    
    % 计算广义特征值
    % B = A0_hat^-1 * A1_hat
    % 检查 A0 是否正定
    try
        chol(A0_hat); 
    catch
        warning('A0 is NOT Positive Definite!');
    end
    
    System_Mat = A0_hat \ A1_hat; % A0^-1 * A1
    eigs_val = eig(System_Mat);
    
    figure;
    plot(real(eigs_val), imag(eigs_val), 'bo', 'MarkerFaceColor', 'b');
    yline(0, 'k--');
    title('Eigenvalues of gPC-SG System');
    xlabel('Real Part'); ylabel('Imaginary Part');
    grid on;
    % exportgraphics(gcf, 'report/hyper4_1.png', 'Resolution', 300);
    % 检查虚部最大值
    max_imag = max(abs(imag(eigs_val)));
    fprintf('Max Imaginary Part: %e\n', max_imag);
    if max_imag < 1e-10
        fprintf('Hyperbolicity Verified: All eigenvalues are real.\n');
    else
        fprintf('Hyperbolicity VIOLATED!\n');
    end
end