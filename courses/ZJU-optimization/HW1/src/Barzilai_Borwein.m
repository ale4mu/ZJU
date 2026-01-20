function [x, val, history] = Barzilai_Borwein(func, grad, x0,parameters)
    % STEEPEST_DESCENT
    %   Inputs:
    %       func         - 目标函数
    %       grad         - 目标函数的梯度
    %       x0           - 初始点
    %       parameters   - 参数
    %   Outputs:
    %       x         - 最优解
    %       val       - 最小值
    %       history   - 迭代历史
    
        % 设置默认参数
        if ~isfield(parameters, 'tol'), parameters.tol = 1e-6; end
        if ~isfield(parameters, 'max_iter'), parameters.max_iter = 1000; end
    
        x = x0;
        history.x_path = {x};
        history.f_path = {func(x)};
        history.alphas = [];
        x_prev = x;
        g_prev = grad(x);
        for i = 1:parameters.max_iter
            g = grad(x);
            
            p = -g;
            if norm(g) < parameters.tol
                fprintf('\nConverged successfully.\n');
                break;
            end
            
            if(i == 1)
                % 回溯线搜索(满足Armigo condition)
                alpha = 1;
                beta = 0.5;
                c = 1e-4;
                val = func(x);
                while func(x + alpha * p) > val + c * alpha * (g' * p)  % Armijo condition
                    alpha = beta * alpha;
                end
            else
                s = x - x_prev;
                y = g - g_prev;
                alpha = (s' * s) / (s' * y);
            end

            x_prev = x;
            g_prev = g;

            % 更新参数
            x = x + alpha * p;
            history.x_path{end+1} = x;
            history.f_path{end+1} = func(x);
            history.alphas(end+1) = alpha;
        end
    
        if i == parameters.max_iter
            fprintf('\nReached max iterations without convergence.\n');
        end
        
        val = func(x);
        history.iterations = i;
    end