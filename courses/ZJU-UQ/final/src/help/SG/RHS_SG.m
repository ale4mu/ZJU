function RHS = RHS_SG(U_hat, dx, gamma, Phi_val, nodes,weights, norm_sq,current_time)

[nVar, nModes, Nx] = size(U_hat); % nModes：M+1, Nx：随机空间网格数
N_quad = length(weights); % 高斯积分点个数

% U_L(:, :, i) 是第 i-1 个单元的右边界值
% U_R(:, :, i) 是第 i 个单元的左边界值
%[U_L, U_R] = WENO_gPC(U_hat, Phi_val, gamma); %example 4.1
%[U_L, U_R] = WENO_gPC43(U_hat, Phi_val, gamma,current_time,nodes,weights,norm_sq); % example 4.3
[U_L, U_R] = WENO_gPC_2D(U_hat, Phi_val, gamma,nodes,weights,norm_sq); % example 4.5

B_plus  = zeros(nVar, nModes, Nx+1); % B正(U+ - U-)
B_minus = zeros(nVar, nModes, Nx+1); % B负(U+ - U-)

% 利用积分路径上的高斯求积：U_L + s * (U_R - U_L)
sqrt_06 = sqrt(0.6);
s_nodes = [-sqrt_06, 0, sqrt_06];
s_weights = [5/9, 8/9, 5/9];
% 映射到[0,1]区间
path_nodes = 0.5 * (s_nodes + 1);
path_weights = 0.5 * s_weights;

% 预计算投影矩阵 P_Mat(q, k) = weights(q) * Phi_val(k, q) / norm_sq(k)
% 用于从物理空间投影回去
P_mat = (weights .* Phi_val') ./ norm_sq'; % (N_quad x nModes)

% % 预计算 Kronecker 积需要的权重矩阵w_q * Phi_q * Phi_q^T
% Kp_cell = cell(N_quad, 1);
% for q = 1:N_quad
%     vec_phi = Phi_val(:, q);
% end
% 下面直接使用了伪谱法，在物理空间乘积再做投影

parfor i = 1:Nx+1
    % i对应论文中的j+1/2
    u_l_coeffs = U_L(:, :, i); % U_{j+1/2}^-
    u_r_coeffs = U_R(:, :, i); % U_{j+1/2}^+
    jump_coeffs = u_r_coeffs - u_l_coeffs; % (U^+ - U^-)

    % 先在物理空间上计算，后面投影回去
    u_l_phys = u_l_coeffs * Phi_val; % (3, N_quad)
    u_r_phys = u_r_coeffs * Phi_val;
    jump_phys = u_r_phys - u_l_phys; %(U_R(xi) - U_L(xi))
    jump_phys_3D = reshape(jump_phys, nVar, 1, N_quad);

    % B = A0^-1 * A1
    % 通过B(path) * Jump(path) = A0(path)^-1 * A1(path) * Jump(path)来计算
    % 下面计算 <Phi_k, B(path) * Jump(path)>
    B_Jump_coeffs = zeros(nVar, nModes);
    local_max_alpha = 0;

    % 用3个gauss积分节点
    for p = 1:3
        s = path_nodes(p);
        w_path = path_weights(p);

        % 三点积分：U(s, xi) = U_L(xi) + s * (U_R(xi) - U_L(xi))
        u_path_phys = u_l_phys + s * jump_phys;

        % Term_phys = zeros(3, N_quad);
        % 下面是循环版本的写法
        % for q = 1:N_quad
        %     u_val = u_path_phys(:, q);
        %     [A_jac, ~, local_alpha] = Euler_Eigen(u_val, gamma); % 求jacobian矩阵和特征值
        %     local_max_alpha = max(local_max_alpha, local_alpha); % 找最大特征值
        %     % B * Jump = A(u) * Jump再做投影
        %     % 对于欧拉方程，B 矩阵 =  对于随机样本xi，计算Jacobian 矩阵 A(U(xi))再投影回系数空间
        %     Term_phys(:, q) = A_jac * jump_phys(:, q);
        % end

        % 对于欧拉方程，B 矩阵 =  对于随机样本xi，计算Jacobian 矩阵 A(U(xi))再投影回系数空间
        % 批量计算所有gauss积分点的jacobian矩阵αF/αU，加快运行速度
        %[A_batch, alpha_vec] = Euler_Eigen_Batch(u_path_phys, gamma); % Example 4.1 4.3
        [A_batch, alpha_vec] = Euler_Eigen_Batch_2D(u_path_phys, gamma); % Example 4.5
        local_max_alpha = max(local_max_alpha, max(alpha_vec)); % 找最大特征值

        % 批量矩阵乘法 Term_phys = A(u) * jump
        % A_batch: (3, 3, N_quad), jump_phys: (3, N_quad) -> reshape to (3, 1, N_quad)
        % Term_phys(:, q) = A_batch(:, :, q) * jump_phys(:, q)
        Term_phys_3D = pagemtimes(A_batch, jump_phys_3D);
        Term_phys = reshape(Term_phys_3D, nVar, N_quad);

        % 投影回系数空间并累加到路径积分
        % Term_phys (3, N_quad) * P_mat (N_quad, nModes) -> (3, nModes)
        Projected = Term_phys * P_mat; % 矩阵乘法版本的投影

        % 注意 weights 是列向量，Phi_val 是 (M+1, N_quad)
        %Projected = (Term_phys * (weights .* Phi_val')) ./ norm_sq';  % 循环版本的投影

        B_Jump_coeffs = B_Jump_coeffs + w_path * Projected;
    end

    % 计算界面通量B正(U^+ - U^-)和B负(U^+ - U^-)
    % 即Lax-Friedrichs 分裂
    B_plus(:, :, i)  = 0.5 * (B_Jump_coeffs + local_max_alpha * jump_coeffs);
    B_minus(:, :, i) = 0.5 * (B_Jump_coeffs - local_max_alpha * jump_coeffs);
end

% 计算RHS的第二部分
Volume_Term_All = zeros(nVar, nModes, Nx);

% P2重构
% P(xi) = a + b*xi + c*xi^2
% 约束:
% 1. int_{-1}^1 P(xi) dxi = 2 * U_avg
% 2. P(-1) = U_left_bound
% 3. P(1)  = U_right_bound
% xi_GL = [-1, 0, 1];
% w_GL  = [1/6, 4/6, 1/6];

parfor j = 1:Nx
    % 论文公式 (2.8) 第二行： - sum w_m B(U) dU/dx
    % 严格来说应该用高斯点处的 WENO 重构值
    % 近似计算：使用单元平均值计算 B，使用重构边界值计算 dU/dx

    % U_R(:, :, j) 是单元 j 左侧重构值 (U_{j-1/2}^+)
    % U_L(:, :, j+1) 是单元 j 右侧重构值 (U_{j+1/2}^-)
    U_left_bound = U_R(:, :, j);
    U_right_bound = U_L(:, :, j+1);
    % U_avg = U_hat(:, :, j);

    dU_dx_coeffs = (U_right_bound - U_left_bound) / dx; % du/dx
    U_center_coeffs = 0.5 * (U_left_bound + U_right_bound);% 这里直接取平均值了
    % 还是投影到物理空间
    u_center_phys = U_center_coeffs * Phi_val;
    du_dx_phys    = dU_dx_coeffs * Phi_val;
    du_dx_phys_3D = reshape(du_dx_phys, nVar, 1, N_quad);

    %Volume_Term_coeffs = zeros(3, nModes);
    % for m = 1:3
    % xi = xi_GL(m);
    % wm = w_GL(m);

    % diff_1 = U_right_bound - U_left_bound;
    % diff_2 = U_right_bound - 2*U_avg + U_left_bound;
    % val_coeffs = U_avg + 0.5*diff_1*xi + 1.5*diff_2*(xi^2 - 1/3);
    % deriv_coeffs = (2/dx) * (0.5*diff_1 + 3*diff_2*xi);
    % val_phys = val_coeffs * Phi_val;
    % deriv_phys = deriv_coeffs * Phi_val;
    % deriv_phys_3D = reshape(deriv_phys, 3, 1, N_quad);

    % 循环版本的写法
    % for q = 1:N_quad
    %     u_val = u_center_phys(:, q);
    %     du_val = du_dx_phys(:, q);

    %     [A_jac, ~, ~] = Euler_Eigen(u_val, gamma);
    %     % 计算B * dU/dx
    %     % 欧拉方程中B = Jocabian矩阵A，因此计算A(u) * dU/dx
    %     term_phys = A_jac * du_val;

    %     % 投影回去
    %     Projected = (term_phys * (weights(q) * Phi_val(:, q)')) ./ norm_sq';
    %     Volume_Term_coeffs = Volume_Term_coeffs + Projected;
    % end

    % 矩阵乘法版本
    % [A_batch, ~] = Euler_Eigen_Batch(val_phys, gamma); % P2重构版本
    % term_phys_3D = pagemtimes(A_batch, deriv_phys_3D);

    %[A_batch, ~] = Euler_Eigen_Batch(u_center_phys, gamma); % 直接用平均值的简化版本
    [A_batch, ~] = Euler_Eigen_Batch_2D(u_center_phys, gamma); % Example 4.5

    term_phys_3D = pagemtimes(A_batch, du_dx_phys_3D);
    term_phys = reshape(term_phys_3D, nVar, N_quad);
    Projected = term_phys * P_mat;
    %Volume_Term_coeffs = Volume_Term_coeffs + wm * Projected;
    Volume_Term_coeffs = Projected;
    %end
    Volume_Term_All(:, :, j) = Volume_Term_coeffs;
end

% 论文公式2.8
B_term = (B_plus(:, :, 1:Nx) +  B_minus(:, :, 2:Nx+1)) / dx; % 前半部分
RHS = -B_term - Volume_Term_All; % 减去后半部分
end











