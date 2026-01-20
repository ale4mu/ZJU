function [Mean_Ref, Std_Ref] = collocaction45(cfg)

    num_sc_points = 30; 
    [nodes_sc, weights_sc] = lgwt(num_sc_points, -1, 1);
    
    Sum_Mean = zeros(cfg.Nx, cfg.Ny);
    Sum_Var  = zeros(cfg.Nx, cfg.Ny);
    
    parfor q = 1:num_sc_points
        xi = nodes_sc(q);
        w_q = weights_sc(q) * 0.5;       
        %rho_final = solve_2D_45(cfg, xi);%example 4.5
        rho_final = solve_2D_46(cfg, xi);
        
        Sum_Mean = Sum_Mean + rho_final * w_q;
        Sum_Var  = Sum_Var  + (rho_final.^2) * w_q;
    end
    
    Mean_Ref = Sum_Mean;
    Var_Diff = Sum_Var - Sum_Mean.^2; 
    Var_Diff(Var_Diff < 0) = 0; 
    Std_Ref = sqrt(Var_Diff);
end