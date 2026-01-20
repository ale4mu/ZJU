function myplot(func,history, x0, x_star)
    % 可视化
    %
    % 输入参数：
    %   history    - 迭代历史
    %   x0         - 初始点
    %   x_star     - 最终最优解
    %   func:      - 目标函数 

    x1_range = linspace(-1.5, 1.5, 200);
    x2_range = linspace(-0.5, 1.5, 200);
    [X1, X2] = meshgrid(x1_range, x2_range);

    Z = arrayfun(@(x1, x2) func([x1; x2]), X1, X2);

    contour(X1, X2, Z, logspace(0, 3.5, 20));
    hold on; grid on; colorbar;

    path = cell2mat(history.x_path);
    plot(path(1, :), path(2, :), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    plot(x0(1), x0(2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Start Point');
    plot(x_star(1), x_star(2), 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'End Point');

    xlabel('x_1'); ylabel('x_2');
    legend('Contours', 'Optimization Path', 'Start Point', 'End Point', 'Location', 'northwest');
    hold off;
end