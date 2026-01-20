function U_new = Euler_forward(U, dt, dx, gamma)
    % 1阶时间离散
    RHS = upwind_SC45(U, dx, gamma);
    U_new = U + dt * RHS;
end