function finiteSystemGroundState(Ls, filename, us, vs, alphas)
    chargeRange = -5:5;
%     s = zeros(length(chargeRange), length(Ls));
%     p = zeros(length(chargeRange), length(Ls));
%     sigma2 = zeros(1, length(Ls));
%     sFull = zeros(1, length(Ls)); 
    ratios = vs - us.';
    ratios(ratios < 0) = 0;
    s = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    p = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    sAlpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s2 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s3 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s4 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s5 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    sigma2 = zeros(length(Ls), length(us), length(vs));
    sFull = zeros(length(Ls), length(us), length(vs)); 
    gaussian = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    for l = 1:length(Ls)
        for l1 = 1:length(us)
            for l2 = 1:length(vs)
                u = us(l1);
                v = vs(l2);
                if (v > u)
                    L = Ls(l);
                    x = L * (v - u) / 2 + chargeRange;
                    cicj = getCiCj0Matrix(L);
                    % Here I trust myself to assign proper u and v values such
                    % that L*u and L*v are integers! This assignment is for
                    % numerical reasons.
                    [~, eigenVals] = eig(cicj(int16(L*u+1):int16(L*v), int16(L*u+1):int16(L*v)));
                    f = zeros(1, length(eigenVals));
                    for i = 1 : length(eigenVals)
                        f(i) = eigenVals(i, i);
                    end
                    p(:, l, l1, l2)  = getSNA(1, f, x, L);
                    s(:, l, l1, l2)  = getEE(f, x, L);
%                     s2(:, l, l1, l2)  = getSNA(2, f, x, L);
%                     s3(:, l, l1, l2)  = getSNA(3, f, x, L);
%                     s4(:, l, l1, l2)  = getSNA(4, f, x, L);
%                     s5(:, l, l1, l2)  = getSNA(5, f, x, L);
                    sAlpha(:, l, l1, l2)  = getSAlpha(1, alphas, f);
                    sFull(l, l1, l2) = sum(s(:, l, l1, l2));
                    fg = fit(x.', real(p(:, l, l1, l2)), gaussian, 'StartPoint', [L * (v - u) / 2 0.1]);
                    sigma2(l, l1, l2) = fg.c;
                end
            end
        end
%     end
%     save(filename, 'Ls', 'p', 's', 'sFull', 'sigma2');
    save(filename, 'p', 's', 's2', 's3', 's4', 's5', 'sAlpha', 'sFull', 'sigma2', 'ratios', 'us', 'vs', 'Ls', 'alphas');
end
