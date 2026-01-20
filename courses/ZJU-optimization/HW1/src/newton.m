function [x,val, history] = newton(func, grad, hessian, x0, parameters)
    % NEWTON法
    %   Inputs:
    %       func         - 目标函数
    %       grad         - 目标函数的梯度
    %       hessian      - 目标函数的Hessian矩阵
    %       x0           - 初始点
    %       parameters   - 参数
    %   Outputs:
    %       x         - 最优解
    %       val       - 最小值
    %       history   - 迭代历史
    
        % 设置默认参数
        if ~isfield(parameters, 'tol'), parameters.tol = 1e-6; end
        if ~isfield(parameters, 'max_iter'), parameters.max_iter = 100; end
    
        x = x0;
        % history.x_path = {x};
        n = length(x0);
        history.x_path = zeros(n, parameters.max_iter + 1); 
        history.x_path(:, 1) = x; 
        history.f_path = {func(x)};
    
        for i = 1:parameters.max_iter
            g = grad(x);
            H = hessian(x);
        
            if norm(g) < parameters.tol
                fprintf('\nConverged successfully.\n');
                break;
            end
    
            % 求解步长 p: H*p = -g
            p = -H \ g;
    
            x = x + p;
            % history.x_path{end+1} = x;
            history.x_path(:, i+1) = x;
            history.f_path{end+1} = func(x);
        end
    
        if i == parameters.max_iter
            fprintf('\nReached max iterations without convergence.\n');
        end
        
        val = func(x);
        history.iterations = i;
        history.x_path = history.x_path(:, 1:i+1);
    end