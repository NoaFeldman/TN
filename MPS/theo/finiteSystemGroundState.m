function finiteSystemGroundState(Ls, filename, ratio1s, ratio2s)
    chargeRange = -5:5;
%     s = zeros(length(chargeRange), length(Ls));
%     p = zeros(length(chargeRange), length(Ls));
%     sigma2 = zeros(1, length(Ls));
%     sFull = zeros(1, length(Ls)); 
    ratios = ratio2s - ratio1s.';
    ratios(ratios < 0) = 0;
    s = zeros(length(chargeRange), length(ratio1s), length(ratio2s));
    p = zeros(length(chargeRange), length(ratio1s), length(ratio2s));
    sigma2 = zeros(length(ratio1s), length(ratio2s));
    sFull = zeros(length(ratio1s), length(ratio2s)); 
    gaussian = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
%     for l = 1:length(Ls)
        for l1 = 1:length(ratio1s)
            for l2 = 1:length(ratio2s)
                ratio1 = ratio1s(l1);
                ratio2 = ratio2s(l2);
                if (ratio2 > ratio1)
        %             L = Ls(l);
                    L = Ls(1);
                    x = L * (ratio2 - ratio1) / 2 + chargeRange;
                    cicj = getCiCj0Matrix(L);
                    [~, v] = eig(cicj((L*ratio1+1):L*ratio2, (L*ratio1+1):L*ratio2));
                    f = zeros(1, length(v));
                    for i = 1 : length(v)
                        f(i) = v(i, i);
                    end
                    p(:, l1, l2)  = getSNA(1, f, x, L);
                    s(:, l1, l2)  = getEE(f, x, L);
                    sFull(l1, l2) = sum(s(:, l1, l2));
                    fg = fit(x.', real(p(:, l1, l2)), gaussian, 'StartPoint', [L * (ratio2 - ratio1) / 2 0.1]);
                    sigma2(l1, l2) = fg.c;
                end
            end
        end
%     end
%     save(filename, 'Ls', 'p', 's', 'sFull', 'sigma2');
    save(filename, 'p', 's', 'sFull', 'sigma2', 'ratios', 'ratio1s', 'ratio2s');
end
