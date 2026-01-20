function order = order(error1, error2, N1, N2)

if error1 == 0
    order = 0.0;
    return;
end

order = log(error1 / error2) / log(N2 / N1);

end