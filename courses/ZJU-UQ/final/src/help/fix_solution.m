function U_hat = fix_solution(U_hat, gamma, Phi_val, epsilon)
    % U_hat: (4, M+1, Nx, Ny)
    [nVar, nModes, Nx, Ny] = size(U_hat);
    
    % 只需要检查均值模态 (k=1)
    % 均值位于 U_hat(:, 1, :, :)
    
    % 提取均值物理量
    rho_mean = squeeze(U_hat(1, 1, :, :));
    m_x_mean = squeeze(U_hat(2, 1, :, :));
    m_y_mean = squeeze(U_hat(3, 1, :, :));
    E_mean   = squeeze(U_hat(4, 1, :, :));
    
    % 计算均值压力
    % 注意：这里用均值动量和均值能量估算压力，这只是一个检查指标
    u_mean = m_x_mean ./ rho_mean;
    v_mean = m_y_mean ./ rho_mean;
    p_mean = (gamma - 1) * (E_mean - 0.5 * rho_mean .* (u_mean.^2 + v_mean.^2));
    
    % 找出坏点
    % 这里的条件要宽松一点，防止误杀，但要足以防止崩溃
    bad_mask = (rho_mean < epsilon) | (p_mean < epsilon);
    
    if any(bad_mask(:))
        fprintf('  Fixing %d bad cells (Mean Violation)...\n', sum(bad_mask(:)));
        
        % 遍历坏点进行修复
        % (由于坏点很少，用循环处理是可以接受的)
        [Ix, Iy] = find(bad_mask);
        
        for k = 1:length(Ix)
            i = Ix(k); j = Iy(k);
            
            % --- Stage 1: 强制修正均值 ---
            % 提取该点的系数向量 (4, M+1)
            u_coeffs = U_hat(:, :, i, j);
            u_0 = u_coeffs(:, 1); % 均值向量
            
            % 1.1 修正密度
            if u_0(1) < epsilon
                u_0(1) = epsilon;
                % 密度太小时，动量清零以防速度爆炸
                u_0(2) = 0; 
                u_0(3) = 0;
            end
            
            % 1.2 修正压力/能量
            rho = u_0(1); mx = u_0(2); my = u_0(3); E = u_0(4);
            p = (gamma-1) * (E - 0.5 * (mx^2 + my^2)/rho);
            
            if p < epsilon
                % 重建能量 E = p_min/(g-1) + 0.5*rho*u^2
                E_new = epsilon/(gamma-1) + 0.5 * (mx^2 + my^2)/rho;
                u_0(4) = E_new;
            end
            
            % --- Stage 2: 抹杀高阶模态 (最稳妥的 Scaling) ---
            % 如果均值都坏了，说明这里的湍流/激波太强，高阶信息已经不可信了。
            % 直接退化为确定性解是防止崩溃的最佳策略。
            % 也就是论文公式中 theta = 0 的情况。
            
            u_coeffs_new = zeros(size(u_coeffs));
            u_coeffs_new(:, 1) = u_0; % 只保留修正后的均值
            
            % 写回 U_hat
            U_hat(:, :, i, j) = u_coeffs_new;
        end
    end
end