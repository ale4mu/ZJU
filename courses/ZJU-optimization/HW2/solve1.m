clear; 
clc; 
close all;

dims = [5, 8, 12, 20]; 
tol = 1e-6;          

figure('Name', 'Q1', 'NumberTitle', 'off');
hold on;
grid on;
set(gca, 'YScale', 'log'); 
title('CG算法在不同维度希尔伯特矩阵上的性能');
xlabel('迭代次数');
ylabel('残差范数');
colors = lines(length(dims)); 

legends = {};

for i = 1:length(dims)
    n = dims(i);
    A = hilb(n); 
    b = ones(n, 1);
    x0 = zeros(n, 1);
    [x, k, residuals] = CG(A, b, x0, tol, 100);
    
    fprintf('维度 n = %d:\n', n);
    fprintf('收敛于 %d 次迭代。\n', k);
    fprintf('最终残差范数: %e\n\n', residuals(end));
    
    iterations = 0:k;
    plot(iterations, residuals, '-o', 'LineWidth', 1.5, 'Color', colors(i,:), 'MarkerSize', 4);
    legends{end+1} = sprintf('n = %d', n);
end

legend(legends, 'Location', 'northeast');
exportgraphics(gcf, "E:\optimization\HW2\Q1.jpg", 'ContentType', 'image');

hold off;