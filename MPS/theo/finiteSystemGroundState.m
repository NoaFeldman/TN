function finiteSystemGroundState(Ls, filename, ratio1, ratio2)
    chargeRange = -5:5;
    s = zeros(length(chargeRange), length(Ls));
    p = zeros(length(chargeRange), length(Ls));
    sigma2 = zeros(1, length(Ls));
    sFull = zeros(1, length(Ls)); 
    gaussian = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    for l = 1:length(Ls)
         L = Ls(l);
         x = L * (ratio2 - ratio1) / 2 + chargeRange;
         cicj = getCiCj0Matrix(L);
         [~, v] = eig(cicj((L*ratio1+1):L*ratio2, (L*ratio1+1):L*ratio2));
         f = zeros(1, length(v));
         for i = 1 : length(v)
             f(i) = v(i, i);
         end
         p(:, l)  = getSNA(1, f, x, L);
         s(:, l)  = getEE(f, x, L);
         sFull(l) = sum(s(:, l));
         fg = fit(x.', real(p(:, l)), gaussian, 'StartPoint', [L * (ratio1 - ratio2) / 2 0.1]);
         sigma2(l) = fg.c;
    end
    save(filename, 'Ls', 'p', 's', 'sFull', 'sigma2');
end
