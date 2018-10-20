function dn = scalingDimension(n, alpha, K)
    c = 1;
    dn = c/12 * (n - n^(-1)) + K / n * (alpha / (2*pi))^2;
end