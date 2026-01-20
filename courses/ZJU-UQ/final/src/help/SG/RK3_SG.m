function U_new = RK3_SG(U, dt, dx, gamma, Phi_val, w_q,norm_sq,current_time,nodes)
    % U 大小: (nVar, M+1, Nx)
    L1 = RHS_SG(U, dx, gamma,Phi_val, nodes,w_q,norm_sq,current_time);
    U1 = U + dt * L1;
    
    
    L2 = RHS_SG(U1, dx, gamma, Phi_val, nodes,w_q,norm_sq,current_time + dt);
    U2 = 0.75 * U + 0.25 * (U1 + dt * L2);
    
    
    L3 = RHS_SG(U2, dx, gamma, Phi_val,nodes,w_q,norm_sq,current_time + dt * 0.5);
    U_new = (1/3) * U + (2/3) * (U2 + dt * L3);
end
