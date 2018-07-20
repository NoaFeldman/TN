function finiteSystemGroundState(Ls, filename, us, vs)
    chargeRange = -5:5;
%     s = zeros(length(chargeRange), length(Ls));
%     p = zeros(length(chargeRange), length(Ls));
%     sigma2 = zeros(1, length(Ls));
%     sFull = zeros(1, length(Ls)); 
    ratios = vs - us.';
    ratios(ratios < 0) = 0;
    s = zeros(length(chargeRange), length(us), length(vs));
    p = zeros(length(chargeRange), length(us), length(vs));
    sigma2 = zeros(length(us), length(vs));
    sFull = zeros(length(us), length(vs)); 
    gaussian = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
%     for l = 1:length(Ls)
        for l1 = 1:length(us)
            for l2 = 1:length(vs)
                % Here I trust myself to assign proper u and v values such
                % that L*u and L*v are integers! This assignment is for
                % numerical reasons.
                u = us(l1);
                v = vs(l2);
                if (v > u)
        %             L = Ls(l);
                    L = Ls(1);
                    x = L * (v - u) / 2 + chargeRange;
                    cicj = getCiCj0Matrix(L);
                    [~, eigenVals] = eig(cicj(int16(L*u+1):int16(L*v), int16(L*u+1):int16(L*v)));
                    f = zeros(1, length(eigenVals));
                    for i = 1 : length(eigenVals)
                        f(i) = eigenVals(i, i);
                    end
                    p(:, l1, l2)  = getSNA(1, f, x, L);
                    s(:, l1, l2)  = getEE(f, x, L);
                    sFull(l1, l2) = sum(s(:, l1, l2));
                    fg = fit(x.', real(p(:, l1, l2)), gaussian, 'StartPoint', [L * (v - u) / 2 0.1]);
                    sigma2(l1, l2) = fg.c;
                end
            end
        end
%     end
%     save(filename, 'Ls', 'p', 's', 'sFull', 'sigma2');
    save(filename, 'p', 's', 'sFull', 'sigma2', 'ratios', 'us', 'vs');
end
