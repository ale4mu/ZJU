% Scheme 3
function [x,u] = scheme3(N,T,a,k)
x_delta = 2*pi / N; % 空间步长
t_delta = 0.9 * x_delta; % 时间步长
lambda = a * t_delta / x_delta;

if a > 0
    % 迭代矩阵
    A = eye(N+1) + lambda * (diag(ones(N,1), 1) - eye(N+1));
    A(N+1,2) = lambda; % 周期边界

    x = linspace(0,2*pi,N+1); % 空间节点
    u = initial(x,k); % 初值条件
    n = floor(T / t_delta); % 时间步数

    if n > 0
        for i = 1:n
            u = u * A'; % Scheme2
        end
    end

    lambda2 = a * (T - n*t_delta) / x_delta;
    % 迭代矩阵2
    A2 = eye(N+1) + lambda2 * (diag(ones(N,1), 1) - eye(N+1));
    A2(N+1,2) = lambda2; % 周期边界

    u = u*A2';
else if a < 0
        A = eye(N+1) + lambda * (eye(N+1) - diag(ones(N,1), -1));
        A(1,N) = -lambda; % 周期边界

        x = linspace(0,2*pi,N+1); % 空间节点
        u = initial(x,k); % 初值条件
        n = floor(T / t_delta); % 时间步数

        if n > 0
            for i = 1:n
                u = u * A'; % Scheme2
            end
        end

        lambda2 = a * (T - n*t_delta) / x_delta;
        % 迭代矩阵2
        A2 = eye(N+1) + lambda2 * (eye(N+1) - diag(ones(N,1), -1));
        A2(1,N) = -lambda2; % 周期边界

        u = u*A2';
end

end