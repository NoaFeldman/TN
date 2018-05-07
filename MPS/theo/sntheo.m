function sn = sntheo(n, t, NA, epsP, epsS, const)
    cn = 1;
    c = 1;
    K = 1;
    x1 = (epsP^2 + t.^2) ./ (epsP / 2);
    x2 = (epsS^2 + t.^2) ./ (epsS / 2);
   
    sn = cn .* x1.^(-c * (n - 1/n) / 6) ...
        .* sqrt(pi * n ./ (2 * K .*  log(x2))) ...
        .* exp(-pi^2 * n * NA^2 ./ (2 * K .* log(x2)));
end