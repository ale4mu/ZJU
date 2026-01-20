function U_out = f_direction(U_in, dt, direction, dx, gamma, Phi_val, weights, norm_sq, nodes,use_low_order)
    [nVar, nModes, Nx, Ny] = size(U_in);
    U_out = U_in;
    
    if strcmp(direction, 'X')
        % X 方向：对每一行 y_j 进行 1D 演化
        parfor j = 1:Ny
            U_slice_1 = U_in(:, :, :, j);
            U_slice = reshape(U_slice_1, nVar, nModes, Nx);
            if use_low_order
                % 调用 1阶求解器 (输入必须是 4xNx 矩阵，去掉中间的 1 维)
                U_slice_2D = squeeze(U_slice); 
                U_new_2D = Euler_forward(U_slice_2D, dt, dx, gamma);
                U_slice_new = reshape(U_new_2D, nVar, 1, Nx);
            else
                U_slice_new = RK3_SG(U_slice, dt, dx, gamma, Phi_val, weights, norm_sq, 0, nodes);
            end
            U_out(:, :, :, j) = U_slice_new;
        end
        
    elseif strcmp(direction, 'Y')
        % Y 方向：对每一列 x_i 进行 1D 演化
        % 原始变量: U = [rho, rho*u, rho*v, E]
        % 交换 u 和 v (第2和第3分量)
        % U' = [rho, rho*v, rho*u, E]
        % 对 U' 求 X 方向通量 F(U') = [rho*v, rho*v^2+p, rho*v*u, v(E+p)] 
        parfor i = 1:Nx
                U_slice_1 = U_in(:, :, i, :); % (4, M+1, Ny)
                U_slice = reshape(U_slice_1, nVar, nModes, Ny);
                U_temp = U_slice;
                U_temp(2, :, :) = U_slice(3, :, :);
                U_temp(3, :, :) = U_slice(2, :, :); 
            if use_low_order
                U_temp_2D = reshape(U_temp, nVar, Ny); 
                U_new_2D = Euler_forward(U_temp_2D, dt, dx, gamma);
                U_temp_new = reshape(U_new_2D, nVar, nModes, Ny);
                U_slice_new = U_temp_new;
            else
                U_temp_new = RK3_SG(U_temp, dt, dx, gamma, Phi_val, weights, norm_sq, 0, nodes);
                U_slice_new = U_temp_new;
            end
            U_slice_new(2, :, :) = U_temp_new(3, :, :);
            U_slice_new(3, :, :) = U_temp_new(2, :, :);
            U_out(:, :, i, :) = U_slice_new;
        end
    end
end