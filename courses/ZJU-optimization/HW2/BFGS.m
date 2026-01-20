function [x,f_val,k,grad_hist] = BFGS(f,grad,x0,tol,max_iter)
% 输入:
%   f             - 目标函数
%   grad          - 梯度
%   x0            - 初始点
%   tol           - 停止迭代的条件
%   max_iter      - 最大迭代次数
%
% 输出:
%   x             - 计算出的最优解
%   f_val         - 最优解处的函数值
%   k             - 实际执行的迭代次数
%   grad_hist     - 梯度范数历史
    n = length(x0);
    x = x0;
    H = eye(n); 
    
    g = grad(x);
    grad_norm = norm(g);
    grad_hist = zeros(max_iter+1,1);
    grad_hist(1) = grad_norm;

    alpha_init = 1.0; % 初始步长
    c1 = 1e-4;   % Armijo Condition
    rho = 0.5;   % 步长缩减因子
    
    for k = 1:max_iter

        if grad_norm < tol
            break;
        end
        
        % p_k = -H_k * g_k
        p = -H * g;
        
        % 回溯线搜索
        f_x = f(x);
        alpha = alpha_init;
        
        % Armijo Condition
        while f(x + alpha * p) > f_x + c1 * alpha * (g' * p)
            alpha = rho * alpha;
        end
        
        ap = alpha * p;
        x_old = x;
        x = x + ap;
        
        g_old = g;
        g = grad(x);
        y = g - g_old;
        s = x - x_old;
        
        if y' * s > 1e-12 % 避免除以零或负数
            rho_k = 1 / (y' * s);
            I = eye(n);
            term1 = I - rho_k * s * y';
            term2 = I - rho_k * y * s';
            H = ((y' *s) / (y' * y)) * I;
            H = term1 * H * term2 + rho_k * s * s';
        end

        grad_norm = norm(g);
        grad_hist(k+1) = grad_norm;
    end
    
    f_val = f(x);
    
    if k == max_iter && grad_norm >= tol
        warning('BFGS:maxIter', 'BFGS在达到最大迭代次数 %d 后仍未收敛。', max_iter);
    end
end