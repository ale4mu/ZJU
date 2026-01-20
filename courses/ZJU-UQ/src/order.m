%% 6. 计算收敛阶
function ord = order(error_prev, error_curr, N_prev, N_curr)
    % 计算收敛阶
    if error_prev == 0 || N_prev == 0
        ord = 0;
    else
        ord = log(error_prev / error_curr) / log(N_curr / N_prev);
    end
end