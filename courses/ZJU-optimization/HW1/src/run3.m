% run3.m

clc;
clear;
close all;

folderpath = 'E:\optimization\HW1';

% Rosenbrock 函数及其导数
rosenbrock = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
grad_rosenbrock = @(x) [
    -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
    200*(x(2) - x(1)^2)
    ];


hess_rosenbrock = @(x) [
    1200*x(1)^2 - 400*x(2) + 2,   -400*x(1);
    -400*x(1),                    200
    ];

% Beale 函数及其导数
beale = @(x) (1.5 - x(1) + x(1)*x(2))^2 + ...
    (2.25 - x(1) + x(1)*x(2)^2)^2 + ...
    (2.625 - x(1) + x(1)*x(2)^3)^2;

grad_beale = @(x) [
    2*((1.5 - x(1) + x(1)*x(2))*(-1+x(2)) + (2.25 - x(1) + x(1)*x(2)^2)*(-1+x(2)^2) + (2.625 - x(1) + x(1)*x(2)^3)*(-1+x(2)^3));
    2*((1.5 - x(1) + x(1)*x(2))*x(1) + (2.25 - x(1) + x(1)*x(2)^2)*2*x(1)*x(2) + (2.625 - x(1) + x(1)*x(2)^3)*3*x(1)*x(2)^2)
    ];

hess_beale = @(x) 2 * [
    (-1+x(2))^2 + (-1+x(2)^2)^2 + (-1+x(2)^3)^2, ...
    (1.5-x(1)*(1-x(2))) + x(1)*(-1+x(2)) + 2*x(2)*(2.25-x(1)*(1-x(2)^2)) + 2*x(1)*x(2)*(-1+x(2)^2) + 3*x(2)^2*(2.625-x(1)*(1-x(2)^3)) + 3*x(1)*x(2)^2*(-1+x(2)^3);

    (1.5-x(1)*(1-x(2))) + x(1)*(-1+x(2)) + 2*x(2)*(2.25-x(1)*(1-x(2)^2)) + 2*x(1)*x(2)*(-1+x(2)^2) + 3*x(2)^2*(2.625-x(1)*(1-x(2)^3)) + 3*x(1)*x(2)^2*(-1+x(2)^3), ...
    x(1)^2 + 2*x(1)*(2.25-x(1)*(1-x(2)^2)) + (2*x(1)*x(2))^2 + 6*x(1)*x(2)*(2.625-x(1)*(1-x(2)^3)) + (3*x(1)*x(2)^2)^2
    ];

%  用最速下降法求Rosenbrock函数的最优解
fprintf('Rosenbrock(steepest descent)\n');
x0 = [-1.2; 1];
parameters.tol = 1e-6;
parameters.max_iter = 100;

[x_star, f_star, history] = steepest_descent(rosenbrock, grad_rosenbrock, x0, parameters);

fprintf('\n--- 结果 ---\n');
fprintf('初始点: x0 = [%.2f, %.2f]\n', x0(1), x0(2));
fprintf('迭代次数: %d\n', history.iterations);
fprintf('最优解: x* = [%.6f, %.6f]\n', x_star(1), x_star(2));
fprintf('最小值: f(x*) = %.6f\n', f_star);
fprintf('最终梯度范数: ||∇f(x*)|| = %.2e\n\n', norm(grad_rosenbrock(x_star)));

figure;
myplot(rosenbrock,history,x0,x_star);
exportgraphics(gcf, fullfile(folderpath, 'rosenbrock_steepest_descent.jpg'), 'ContentType', 'image');

%  用牛顿法求Rosenbrock函数的最优解
fprintf('Rosenbrock(newton)\n');

[x_star, f_star, history] = newton(rosenbrock, grad_rosenbrock,hess_rosenbrock, x0, parameters);

fprintf('\n--- 结果 ---\n');
fprintf('初始点: x0 = [%.2f, %.2f]\n', x0(1), x0(2));
fprintf('迭代次数: %d\n', history.iterations);
fprintf('最优解: x* = [%.6f, %.6f]\n', x_star(1), x_star(2));
fprintf('最小值: f(x*) = %.6f\n', f_star);
fprintf('最终梯度范数: ||∇f(x*)|| = %.2e\n\n', norm(grad_rosenbrock(x_star)));

figure;
myplot(rosenbrock,history,x0,x_star);
exportgraphics(gcf, fullfile(folderpath, 'rosenbrock_newton.jpg'), 'ContentType', 'image');

%  用最速下降法求Beale函数的最优解
fprintf('Beale(steepest descent)\n');
x0 = [-1.2; 1]; 
parameters_nm.tol = 1e-6;
parameters_nm.max_iter = 100;

[x_star, f_star, history] = steepest_descent(beale, grad_beale, x0, parameters_nm);

fprintf('\n--- 结果 ---\n');
fprintf('初始点: x0 = [%.2f, %.2f]\n', x0(1), x0(2));
fprintf('迭代次数: %d\n', history.iterations);
fprintf('最优解: x* = [%.6f, %.6f]\n', x_star(1), x_star(2));
fprintf('最小值: f(x*) = %.6f\n', f_star);
fprintf('最终梯度范数: ||∇f(x*)|| = %.2e\n\n', norm(grad_beale(x_star)));

figure;
myplot(beale,history,x0,x_star);
exportgraphics(gcf, fullfile(folderpath, 'beale_steepest_descent.jpg'), 'ContentType', 'image');



%  用牛顿法求Beale函数的最优解
fprintf('Beale(newton)\n');
x0 = [-1.2; 1]; 
parameters_nm.tol = 1e-6;
parameters_nm.max_iter = 100;

[x_star, f_star, history] = newton(beale, grad_beale, hess_beale, x0, parameters_nm);

fprintf('\n--- 结果 ---\n');
fprintf('初始点: x0 = [%.2f, %.2f]\n', x0(1), x0(2));
fprintf('迭代次数: %d\n', history.iterations);
fprintf('最优解: x* = [%.6f, %.6f]\n', x_star(1), x_star(2));
fprintf('最小值: f(x*) = %.6f\n', f_star);
fprintf('最终梯度范数: ||∇f(x*)|| = %.2e\n\n', norm(grad_beale(x_star)));

figure;
myplot(beale,history,x0,x_star);
exportgraphics(gcf, fullfile(folderpath,'beale_newton.jpg'), 'ContentType', 'image');


