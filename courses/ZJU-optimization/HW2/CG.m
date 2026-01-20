function [x,k,residuals] = CG(A,b,x0,tol,max_iter)

% 输出:
%   x          - 解
%   k          - 迭代次数
%   residuals  - 残差
    x = x0;
    r = A * x0 -b;
    p = -r;
    residual = norm(r);
    residuals = zeros(max_iter+1, 1);
    residuals(1) = residual;

    if residual < tol
        k = 0;
        return;
    end

    for k = 1:max_iter
        Ap = A * p;
        alpha =  - (r' * p) / (p' * Ap);
        x = x + alpha * p;
        r_old = r;
        r = r + alpha * Ap;
        residual = norm(r);
        residuals(k+1) = residual;
        if residual < tol
            break;
        end
        beta = (r' * r) / (r_old' * r_old);
        p = -r + beta * p;
    end

    if k == max_iter && residual >= tol
        warning('CG:maxIter', '共轭梯度算法在达到最大迭代次数 %d 后仍未收敛。', max_iter);
    end

    residuals = residuals(1:k+1);
end