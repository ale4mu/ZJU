% 检查是否超出容许状态集
function [is_valid, min_rho,min_p] = check_Euler(U_coeffs, Phi_val, gamma)
    % U_coeffs: (3, nModes)
    % Phi_val: (nModes, N_quad)

    is_valid = true;
    U_phys = U_coeffs * Phi_val;
    rho = U_phys(1, :);
    m   = U_phys(2, :);
    E   = U_phys(3, :);
    % 检查密度
    min_rho = min(rho);
    % 检查压力
    u = m ./ rho;
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);
    min_p = min(p);
    if min_p <= 1e-2 || min_rho <= 1e-2
        is_valid = false;
    end
end