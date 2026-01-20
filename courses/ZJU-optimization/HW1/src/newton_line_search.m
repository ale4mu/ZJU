function [x, val, history] = newton_line_search(func, grad, hessian, x0, parameters)

        if ~isfield(parameters, 'tol'), parameters.tol = 1e-6; end
        if ~isfield(parameters, 'max_iter'), parameters.max_iter = 100; end
    
        x = x0;
        n = length(x0);
        % history.x_path = {x};
        history.x_path = zeros(n, parameters.max_iter + 1); 
        history.x_path(:, 1) = x; 
        history.f_path = {func(x)};
        history.alphas = [];

        for i = 1:parameters.max_iter
            g = grad(x);
            H = hessian(x);
            f_val = func(x);
            
            if norm(g) < parameters.tol
                fprintf('\nConverged successfully.\n');
                break;
            end
            
            % 求解步长 p
            p = -H \ g;
            
            % 回溯线搜索 (满足Armijo Condition)
            alpha = 1;
            % alpha = alpha_init;      % 初始步长
            beta = 0.5;     
            c = 0.6;       
            while func(x + alpha * p) > f_val + c * alpha * (g' * p) % Armijo Condition
                alpha = beta * alpha;
            end
            % if(alpha ~= 1)
            %     fprintf('alpha reduced to %f at iteration %d\n', alpha, i);
            % end
            history.alphas(end+1) = alpha;

            % 更新
            x = x + alpha * p;
            history.x_path(:, i+1) = x;
            % history.x_path{end+1} = x;
            history.f_path{end+1} = func(x);
        end
    
        if i == parameters.max_iter
            fprintf('\nReached max iterations without convergence.\n');
        end
        
        val = func(x);
        history.iterations = i;
        history.x_path = history.x_path(:, 1:i+1);
    end