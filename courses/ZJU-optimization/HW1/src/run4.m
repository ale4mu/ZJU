clc;
clear;
close all;

A = [5   1   0   1/2;
     1   4   1/2 0;
     0   1/2 3   0;
     1/2 0   0   2];

f = @(x, sigma) 0.5*(x'*x) + 0.25*sigma*(x'*A*x)^2;
grad_f = @(x, sigma) x + sigma*(x'*A*x)*(A*x);
hess_f = @(x, sigma) eye(4) + 2*sigma*A*(x*x')*A + sigma*(x'*A*x)*A;

c70 = cosd(70); s70 = sind(70);
x0_1 = [c70; s70; c70; s70];

c50 = cosd(50); s50 = sind(50);
x0_2 = [c50; s50; c50; s50];

sigmas = [1, 1e4];
initial_points = {x0_1, x0_2};
initial_points_names = {'x0=(cos70,sin70,..)', 'x0=(cos50,sin50,..)'};

parameters.tol = 1e-6;
parameters.max_iter = 50;

fprintf('%-8s | %-22s | %-22s | %-10s | %-20s | %-15s\n', ...
    'Sigma', 'Initial Point', 'Method', 'Iterations', 'x*','f(x*)');

for i = 1:length(sigmas)
    sigma = sigmas(i);

    func = @(x) f(x, sigma);
    grad = @(x) grad_f(x, sigma);
    hessian = @(x) hess_f(x, sigma);
    
    for j = 1:length(initial_points)
        x0 = initial_points{j};
        
        [x_pure, f_val_pure, hist_pure] = newton(func, grad, hessian, x0, parameters);
        fprintf('%-8.1e | %-22s | %-22s | %-10d | [%.6f, %.6f, %.6f, %.6f] | %-15.6e\n', ...
            sigma, initial_points_names{j}, 'Pure Newton', ...
            hist_pure.iterations, x_pure(1), x_pure(2), x_pure(3), x_pure(4), f_val_pure);

        % analyze_convergence(hist_pure, 'Pure Newton');

        [x_ls, f_val_ls, hist_ls] = newton_line_search(func, grad, hessian, x0, parameters);
        fprintf('%-8.1e | %-22s | %-22s | %-10d | [%.6f, %.6f, %.6f, %.6f] | %-15.6e\n', ...
            sigma, initial_points_names{j}, 'Newton with Line Search', ...
            hist_ls.iterations,x_ls(1),x_ls(2),x_ls(3),x_ls(4), f_val_ls);

        % analyze_convergence(hist_ls, 'Newton with Line Search');
    end
end