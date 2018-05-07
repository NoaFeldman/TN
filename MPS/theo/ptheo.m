function sn = ptheo(t, NA, epsP, cn)
    n = 1;
    K = 1;
    x = (epsP^2 + t.^2) ./ (epsP / 2);
   
    sn = cn .* sqrt(pi * n ./ (2 * K .*  log(x))) ...
        .* exp(-pi^2 * n * NA^2 ./ (2 * K .* log(x)));
end