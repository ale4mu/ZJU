clc;
clear;
close all;

tol = 1e-5;
max_iter = 500; 
dims = [6, 8, 10];

for n = dims
    fprintf('维度 n = %d\n', n);
    
    x0 = ones(n, 1);
    x0(1:2:n) = -1.2;
    
    tic; 
    [x, f_val, k, grad_hist] = BFGS(@rosenbrock, @rosenbrock_grad, x0, tol, max_iter);
    elapsed_time = toc;
    
    fprintf('BFGS 算法在 %d 次迭代后收敛。\n', k);
    fprintf('耗时: %.4f 秒。\n', elapsed_time);
    fprintf('最终函数值: %e\n', f_val);
    fprintf('最终梯度范数: %e\n', grad_hist(end));
    fprintf('最终解(前4个分量): [%.4f, %.4f, %.4f, %.4f, ...]\n\n', x(1:4));
end

x0 = [3; -1; 0; 1];
n = length(x0);
tic;
[x, f_val, k, grad_hist] = BFGS(@powell, @powell_grad, x0, tol, max_iter);
elapsed_time = toc;

fprintf('BFGS 算法在 %d 次迭代后收敛。\n', k);
fprintf('耗时: %.4f 秒。\n', elapsed_time);
fprintf('最终函数值: %e\n', f_val);
fprintf('最终梯度范数: %e\n', grad_hist(end));
fprintf('最终解: [%.4f, %.4f, %.4f, %.4f]\n\n', x);
