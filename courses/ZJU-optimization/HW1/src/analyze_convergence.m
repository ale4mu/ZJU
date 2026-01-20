function analyze_convergence(history, method_name)
    x_history = history.x_path;
    num_iters = size(x_history, 2) - 1;
    start_iter = max(1, num_iters - 5);
    
    fprintf('--- Convergence Rate Analysis for: %s ---\n', method_name);
    fprintf('Iter (k) |  ||x_k||   | 二次收敛率 | 线性收敛率\n');
    for k = start_iter:num_iters-1
        x_current = x_history(:, k+1); 
        x_next    = x_history(:, k+2); 
        norm_current = norm(x_current);
        norm_next = norm(x_next);
        
        if norm_current < 1e-15
            ratio1 = NaN;
            ratio2 = NaN;
        else
            ratio1 = norm_next / (norm_current^2);
            ratio2 = norm_next / norm_current;
        end
        fprintf('%8d | %.3e | %.3e | %.4f | %.4f\n', k, norm_current, norm_next, ratio1,ratio2);
    end
end