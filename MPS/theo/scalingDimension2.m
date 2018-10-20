function dn2 = scalingDimension2(n, alpha, K)
    if mod(n, 2) == 0
        dn2 = 2 * scalingDimension(n/2, alpha, K);
    else
        dn2 = scalingDimension(n, alpha, K);
    end
end