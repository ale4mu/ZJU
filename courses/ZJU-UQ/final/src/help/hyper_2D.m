function hyper_2D(U_hat_slice, gamma, Phi_val, weights, norm_sq)
    
    [nVar, nModes] = size(U_hat_slice);
    N_quad = length(weights);  

    [is_valid, ~, ~] = check_Euler_2D(U_hat_slice, Phi_val, gamma);
    if ~is_valid
        warning('Invalid solution detected (negative density/pressure)! Hyperbolicity check might fail.');
    end
    
    U_phys = U_hat_slice * Phi_val; % (nVar, N_quad)
    

    A0_hat = zeros(nVar*nModes, nVar*nModes);
    A1_hat = zeros(nVar*nModes, nVar*nModes);
    
    for q = 1:N_quad
        u_val = U_phys(:, q);
        w = weights(q);

        [L, Lambda] = Euler_Eigen_2D(u_val, gamma); 
        
        A0_phys = L' * L;
        A1_phys = L' * Lambda * L;
        
        for k = 1:nModes
            for l = 1:nModes
                r_idx = (k-1)*nVar + (1:nVar);
                c_idx = (l-1)*nVar + (1:nVar);
                
                basis_prod = Phi_val(k, q) * Phi_val(l, q);
                
                A0_hat(r_idx, c_idx) = A0_hat(r_idx, c_idx) + w * A0_phys * basis_prod;
                A1_hat(r_idx, c_idx) = A1_hat(r_idx, c_idx) + w * A1_phys * basis_prod;
            end
        end
    end
    
    try
        chol(A0_hat); 
    catch
        warning('A0 is NOT Positive Definite!');
    end
    
    System_Mat = A0_hat \ A1_hat;
    eigs_val = eig(System_Mat);
    
    figure;
    plot(real(eigs_val), imag(eigs_val), 'bo', 'MarkerFaceColor', 'b');
    yline(0, 'k--');
    title(sprintf('Eigenvalues'));
    xlabel('Real Part'); ylabel('Imaginary Part');
    grid on;
    
    max_imag = max(abs(imag(eigs_val)));
    fprintf('Max Imaginary Part: %e\n', max_imag);
    if max_imag < 1e-8
        fprintf('Hyperbolicity Verified.\n');
    else
        fprintf('Hyperbolicity VIOLATED!\n');
    end
end