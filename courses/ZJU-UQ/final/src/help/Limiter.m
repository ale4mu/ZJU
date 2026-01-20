% 限制器，把超出容许状态集的拉回
function U_limited = Limiter(U_target, U_base, Phi_val, gamma)
    % U_target: 待检查的WENO重构值
    % U_base:   单元平均值
   
    [is_valid,~,~] = check_Euler(U_target, Phi_val, gamma);
    if is_valid
        U_limited = U_target;
        return;
    end
    % 如果超过限制，开始寻找 theta
    % U_new = U_base + theta * (U_target - U_base)
    dU = U_target - U_base;
    % 二分法查找最大 theta
    theta_low = 0;
    theta_high = 1;
    tol = 1e-6; 
    for iter = 1:40
        theta_mid = 0.5 * (theta_low + theta_high);
        U_test = U_base + theta_mid * dU;
        [is_valid,~,~] = check_Euler(U_test, Phi_val, gamma);
        if is_valid
            theta_low = theta_mid; % 可行，尝试更大的 theta
        else
            theta_high = theta_mid; % 不可行，减小 theta
        end
        if (theta_high - theta_low) < tol
            break;
        end
    end
    % 返回修正值
    U_limited = U_base + theta_low * dU;
end
