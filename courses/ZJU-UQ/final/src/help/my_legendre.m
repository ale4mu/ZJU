function val = my_legendre(k,x)
    if k == 0
        val = 1;
    elseif k == 1
        val = x;
    else
        P_prev2 = 1; % P0(x)
        P_prev1 = x; % P1(x)
        for n = 2:k
            P_curr = ((2*n - 1) * x .* P_prev1 - (n - 1) * P_prev2) / n;
            P_prev2 = P_prev1;
            P_prev1 = P_curr;
        end
        val = P_curr;
    end
end