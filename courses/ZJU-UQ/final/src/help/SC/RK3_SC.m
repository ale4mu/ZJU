% Example 4.1
function U_new = RK3_SC(U, dt, dx, gamma)
    % Runge-Kutta
    % U_1 = U_n + dt * L(U_n)
    L_n = WENO(U, dx, gamma); % 全局Lax-Friedrichs 分裂
    %L_n = WENO_LLF(U, dx, gamma); % 局部Lax-Friedrichs 分裂
    U1 = U + dt * L_n;

    % U_2 = 3/4 U_n + 1/4 (U_1 + dt * L(U_1))
    L_1 = WENO(U1, dx, gamma);
    %L_1 = WENO_LLF(U1, dx, gamma);
    U2 = 0.75 * U + 0.25 * (U1 + dt * L_1);

    % U_new = 1/3 U_n + 2/3 (U_2 + dt * L(U_2))
    L_2 = WENO(U2, dx, gamma);
    %L_2 = WENO_LLF(U2, dx, gamma);
    U_new = (1/3) * U + (2/3) * (U2 + dt * L_2);
end
