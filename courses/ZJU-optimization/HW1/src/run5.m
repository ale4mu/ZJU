clc;
clear;
close all;

m = 512;
n = 1024;
r = 0.1; % 稀疏度
mu_values = [1e-2, 1e-3]; % 不同的 mu 值
parameters.tol = 1e-6;
parameters.max_iter = 2000;

for k = 1:length(mu_values)
    mu = mu_values(k);
    delta = 1e-3 * mu;
    
    fprintf('\n======================================================\n');
    fprintf('           Running for mu = %.1e, delta = %.1e\n', mu, delta);
    fprintf('======================================================\n');

    rng(0); 
    A = randn(m, n);
    x_true = sprandn(n, 1, r); 
    b = A * x_true;
    x0 = zeros(n, 1); % 初始点

   
    % 正则项 L_delta(x) 的梯度
    grad_L = @(x) (abs(x) < delta) .* (x / delta) + (abs(x) >= delta) .* sign(x);
    
    % 目标函数的梯度
    grad_f = @(x) A' * (A * x - b) + mu * grad_L(x);
    
    % 目标函数 
    l_delta = @(x) (abs(x) < delta) .* (0.5/delta * x.^2) + (abs(x) >= delta) .* (abs(x) - 0.5*delta);
    L_delta = @(x) sum(l_delta(x));
    func = @(x) 0.5 * norm(A * x - b)^2 + mu * L_delta(x);

    fprintf('\n--- Solving with Steepest Descent Method ---\n');
    tic;
    [x,val,history] = steepest_descent(func, grad_f, x0, parameters);
    time = toc;
    
    err = norm(x - x_true) / norm(x_true);
    
    fprintf('Iterations: %d\n', history.iterations);
    fprintf('Final function value: %.6f\n', history.f_path{end});
    fprintf('||x - x_true|| / ||x_true||: %.6f\n', err);
    fprintf('Time taken: %.4f seconds\n', time);

    fprintf('\n--- Solving with Barzilai-Borwein Method ---\n');
    tic;
    [x_bb,val_bb, history_bb] = Barzilai_Borwein(func,grad_f, x0,parameters);
    time_bb = toc;
    err_bb = norm(x_bb - x_true) / norm(x_true);
    
    fprintf('Iterations: %d\n', history_bb.iterations);
    fprintf('Final function value: %.6f\n', history_bb.f_path{end}); 
    fprintf('||x - x_true|| / ||x_true||: %.6f\n', err_bb);
    fprintf('Time taken: %.4f seconds\n', time_bb);
    
    figure('Name', sprintf('Sparsity Recovery (mu = %.1e)', mu));
    subplot(3, 1, 1);
    stem(x_true, 'k', 'Marker', 'none');
    title('True X');
    % legend('x_{true}');
    grid on;
    
    folderpath = 'E:\optimization\HW1';
    subplot(3, 1, 2);
    stem(x, 'b', 'Marker', 'none');
    title(sprintf('Steepest Descent (Err=%.4f)', err));
    % legend('x_{SD}');
    grid on;
    exportgraphics(gcf, fullfile(folderpath,sprintf('steepest_descent_mu_%.1e.png', mu)), 'ContentType', 'image');


    subplot(3, 1, 3);
    stem(x_bb, 'r', 'Marker', 'none');
    title(sprintf('Barzilai-Borwein (Err=%.4f)', err_bb));
    % legend('x_{BB}');
    grid on;
    exportgraphics(gcf, fullfile(folderpath, sprintf('BB5_mu_%.1e.png', mu)), 'ContentType', 'image');
end