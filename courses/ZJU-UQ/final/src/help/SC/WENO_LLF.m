function RHS = WENO_LLF(U, dx, gamma)
    [nVar, Nx] = size(U);
    
    % 扩展 U 以处理周期性边界 [N-2, N-1, N, 1...N, 1, 2, 3]
    idx_ext = [Nx-2, Nx-1, Nx, 1:Nx, 1, 2, 3];
    U_ext = U(:, idx_ext);
    
    U_L = zeros(nVar, Nx+1); % 左侧值
    U_R = zeros(nVar, Nx+1); % 右侧值
    
    for i = 1:Nx
        idx = i + 3; %前后各加了3个
        
        for k = 1:nVar
            % [i-2, ..., i+2]
            stencil = U_ext(k, idx-2 : idx+2);
            U_L(k, i+1) = WENO_pos(stencil);

            stencil_neg = U_ext(k, idx-1 : idx+3);
            U_R(k, i+1) = WENO_neg(stencil_neg); 
        end
    end
    % 周期性边界处理
    U_L(:, 1) = U_L(:, Nx+1);
    U_R(:, 1) = U_R(:, Nx+1);
    
    F_hat = zeros(nVar, Nx+1);
    
    for i = 1:Nx+1
        Ul = U_L(:, i);
        Ur = U_R(:, i);
        
        % 计算左右状态的物理通量 F(Ul), F(Ur)
        Fl = Flux(Ul, gamma);
        Fr = Flux(Ur, gamma);
        
        % 计算局部最大特征波速 alpha
        al = MaxWaveSpeed(Ul, gamma);
        ar = MaxWaveSpeed(Ur, gamma);
        alpha = max(al, ar);
        
        % LLF 通量公式
        F_hat(:, i) = 0.5 * (Fl + Fr - alpha * (Ur - Ul));
    end
    
    RHS = - (F_hat(:, 2:Nx+1) - F_hat(:, 1:Nx)) / dx;
end

function F = Flux(U, gamma)
    rho = U(1); 
    m = U(2); 
    E = U(3);
    % epsilon = 1e-10;
    u = m / rho;
    p = (gamma - 1) * (E - 0.5 * rho * u^2);
    
    F = [m; m*u + p; (E+p)*u];
end

function a = MaxWaveSpeed(U, gamma)
    rho = U(1); 
    m = U(2); 
    E = U(3);
    % epsilon = 1e-10;
    u = m / rho;
    p = (gamma - 1) * (E - 0.5 * rho * u^2);
    c = sqrt(gamma * p / rho );
    a = abs(u) + c;
end
