% 初始化函数
function U_coeffs = Initial_gPC(x, y, cfg, nodes, weights, Phi_val, norm_sq)
% Example 4.5
% x, y in [0, 1]

% 确定象限
if x > 0.5 && y > 0.5
    rho0 = 1; u0 = 0; v0 = 0; p0 = 1;
elseif x < 0.5 && y > 0.5
    rho0 = 0.5197; u0 = -0.7259; v0 = 0; p0 = 0.4;
elseif x < 0.5 && y < 0.5
    rho0 = 1; u0 = -0.7259; v0 = -0.7259; p0 = 1;
else
    rho0 = 0.5197; u0 = 0; v0 = -0.7259; p0 = 0.4;
end

val = -0.7259;

u_phys = u0 * ones(size(nodes));
if abs(u0 - val) < 1e-6
    u_phys = u0 + 0.1 * nodes;
end

v_phys = v0 * ones(size(nodes));
if abs(v0 - val) < 1e-6
    v_phys = v0 + 0.1 * nodes;
end

rho_phys = rho0 * ones(size(nodes));
p_phys   = p0 * ones(size(nodes));

q2 = u_phys.^2 + v_phys.^2;
E_phys = p_phys/(cfg.gamma-1) + 0.5 * rho_phys .* q2;

% 投影
U_coeffs = zeros(4, cfg.M+1);
for k = 1:(cfg.M+1)
    basis = Phi_val(k, :)';
    % rho
    U_coeffs(1, k) = sum(rho_phys .* weights .* basis) / norm_sq(k);
    % m_x
    U_coeffs(2, k) = sum((rho_phys .* u_phys) .* weights .* basis) / norm_sq(k);
    % m_y
    U_coeffs(3, k) = sum((rho_phys .* v_phys) .* weights .* basis) / norm_sq(k);
    % E
    U_coeffs(4, k) = sum(E_phys .* weights .* basis) / norm_sq(k);
end
end