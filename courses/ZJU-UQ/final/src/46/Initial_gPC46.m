function U_coeffs = Initial_gPC46(x, y, cfg, nodes, weights, Phi_val, norm_sq)
    % Example 4.6 
    % 密度扰动: rho = rho0 * (1 + 0.1*xi)

    if x > 0.5 && y > 0.5
        rho0=0.5197; u0=0.1; v0=0.1; p0=0.4;
    elseif x < 0.5 && y > 0.5
        rho0=1; u0=-0.6259; v0=0.1; p0=1;
    elseif x < 0.5 && y < 0.5
        rho0=0.8; u0=0.1; v0=0.1; p0=1;
    else
        rho0=1; u0=0.1; v0=-0.6259; p0=1;
    end
    
    rho_phys = rho0 * (1 + 0.1 * nodes);

    u_phys = u0 * ones(size(nodes));
    v_phys = v0 * ones(size(nodes));
    p_phys = p0 * ones(size(nodes));
 
    q2 = u_phys.^2 + v_phys.^2;
    E_phys = p_phys / (cfg.gamma - 1) + 0.5 * rho_phys .* q2;
    
    U_coeffs = zeros(4, cfg.M+1);
    for k = 1:(cfg.M+1)
        basis = Phi_val(k, :)';
        % rho
        U_coeffs(1, k) = sum(rho_phys .* weights .* basis) / norm_sq(k);
        % m_x = rho * u
        U_coeffs(2, k) = sum((rho_phys .* u_phys) .* weights .* basis) / norm_sq(k);
        % m_y = rho * v
        U_coeffs(3, k) = sum((rho_phys .* v_phys) .* weights .* basis) / norm_sq(k);
        % E
        U_coeffs(4, k) = sum(E_phys .* weights .* basis) / norm_sq(k);
    end
end