function [is_valid, min_rho, min_p] = check_Euler_2D(U_coeffs, Phi_val, gamma)
    % U_coeffs: (nVar, nModes)
    % Phi_val: (nModes, N_quad)
    
    [nVar, ~] = size(U_coeffs);

    U_phys = U_coeffs * Phi_val; % (nVar, N_quad)
    
    rho = U_phys(1, :);
        
    mx = U_phys(2, :);
    my = U_phys(3, :);
    E  = U_phys(4, :);
    q2 = (mx.^2 + my.^2) ./ (rho.^2); % u^2 + v^2
        
    p = (gamma - 1) * (E - 0.5 * rho .* q2);
    
    min_rho = min(rho);
    min_p   = min(p);
    
    epsilon = 1e-2; 
    if min_rho < epsilon || min_p < epsilon
        is_valid = false;
    else
        is_valid = true;
    end
end