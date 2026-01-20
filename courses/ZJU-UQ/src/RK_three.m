% 三阶Runge−Kutta方法
function [t,u] = RK_three(N,a,b,beta,alpha)
    t_delta = (b-a)/N; % 时间步长

    t = linspace(a,b,N+1); % 节点
    u = zeros(1,length(t)); % 节点函数值

    u(1) = beta; % 初值条件
    for i = 2:length(t)
        % 三阶Runge−Kutta方法
        u_temp1 = u(i-1) - t_delta * alpha * u(i-1);
        u_temp2 = 3/4 * u(i-1) + 1/4 * u_temp1 - 1/4 *t_delta * alpha * u_temp1;
        u(i) = 1/3 * u(i-1) + 2/3 * u_temp2 - 2/3 * t_delta * alpha * u_temp2;
    end
end