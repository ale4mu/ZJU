function val = WENO_pos(v)
    % v 是长度为5的向量: [v1, v2, v3, v4, v5]
    % 对应 [j-2, j-1, j, j+1, j+2]
    % 平滑度因子，beta越小，模板越光滑
    beta1 = 13/12 * (v(1) - 2*v(2) + v(3))^2 + 1/4 * (v(1) - 4*v(2) + 3*v(3))^2;
    beta2 = 13/12 * (v(2) - 2*v(3) + v(4))^2 + 1/4 * (v(2) - v(4))^2;
    beta3 = 13/12 * (v(3) - 2*v(4) + v(5))^2 + 1/4 * (3*v(3) - 4*v(4) + v(5))^2;
    d1 = 0.1; d2 = 0.6; d3 = 0.3;
    eps = 1e-2; % 这里不加保护在example 4.3中会崩溃
    alpha1 = d1 / (beta1+eps)^2; % 非线性权重
    alpha2 = d2 / (beta2+eps)^2;
    alpha3 = d3 / (beta3+eps)^2;
    w_sum = alpha1 + alpha2 + alpha3; % 归一化
    w1 = alpha1 / w_sum;
    w2 = alpha2 / w_sum;
    w3 = alpha3 / w_sum;
    q1 = 1/3*v(1) - 7/6*v(2) + 11/6*v(3);
    q2 = -1/6*v(2) + 5/6*v(3) + 1/3*v(4);
    q3 = 1/3*v(3) + 5/6*v(4) - 1/6*v(5);
    val = w1*q1 + w2*q2 + w3*q3;
end





