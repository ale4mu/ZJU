%% 1. 精确解函数
function [ex_mean, ex_var] = exact_solution(t)
    % 精确解: y(t,α) = β * exp(-α*t)
    % 期望: E[y(t)] = exp(t^2/2)
    % 方差: Var[y(t)] = exp(2t^2) - exp(t^2)
    ex_mean = exp(t.^2 / 2);
    ex_var = exp(2*t.^2) - exp(t.^2);
end