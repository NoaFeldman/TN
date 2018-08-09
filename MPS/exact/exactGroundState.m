function exactGroundState(Ls, filename, us, vs, alphas)
    chargeRange = -5:5;
%     s = zeros(length(chargeRange), length(Ls));
%     p = zeros(length(chargeRange), length(Ls));
%     sigma2 = zeros(1, length(Ls));
%     sFull = zeros(1, length(Ls)); 
    ratios = vs - us.';
    ratios(ratios < 0) = 0;
    s = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s1 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s1Alpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s2Alpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s3Alpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s4Alpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s5Alpha = zeros(length(alphas), length(Ls), length(us), length(vs));
    s2 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s3 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s4 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    s5 = zeros(length(chargeRange), length(Ls), length(us), length(vs));
    varN1 = zeros(length(Ls), length(us), length(vs));
    varN2 = zeros(length(Ls), length(us), length(vs));
    varN3 = zeros(length(Ls), length(us), length(vs));
    varN4 = zeros(length(Ls), length(us), length(vs));
    varN5 = zeros(length(Ls), length(us), length(vs));
    sn02 = zeros(length(Ls), length(us), length(vs));
    sn03 = zeros(length(Ls), length(us), length(vs));
    sn04 = zeros(length(Ls), length(us), length(vs));
    sn05 = zeros(length(Ls), length(us), length(vs));
    sFull = zeros(length(Ls), length(us), length(vs)); 
    gaussian = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    renyi = fittype('a*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
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
                    s(:, l, l1, l2)  = getExactEE(f, x, L);
                    s1(:, l, l1, l2)  = getSNA(1, f, x, L);
                    s2(:, l, l1, l2)  = getSNA(2, f, x, L);
                    s3(:, l, l1, l2)  = getSNA(3, f, x, L);
                    s4(:, l, l1, l2)  = getSNA(4, f, x, L);
                    s5(:, l, l1, l2)  = getSNA(5, f, x, L);
                    s1Alpha(:, l, l1, l2)  = getSAlpha(1, alphas, f, 1);
                    s2Alpha(:, l, l1, l2)  = getSAlpha(2, alphas, f, 1);
                    s3Alpha(:, l, l1, l2)  = getSAlpha(3, alphas, f, 1);
                    s4Alpha(:, l, l1, l2)  = getSAlpha(4, alphas, f, 1);
                    s5Alpha(:, l, l1, l2)  = getSAlpha(5, alphas, f, 1);
                    sFull(l, l1, l2) = sum(s(:, l, l1, l2));

                    sAlpha = getExactEEForFlux(f, alphas, L);
                    sFullAlpha(l2) = sum(sAlpha);
                    
%                     fg = fit(x.', real(s1(:, l, l1, l2)), gaussian, 'StartPoint', [L * (v - u) / 2 0.1]);
%                     varN1(l, l1, l2) = fg.c;
%                     fg = fit(x.', real(s2(:, l, l1, l2)), renyi, 'StartPoint', [1 L * (v - u) / 2 0.1]);
%                     sn02(l, l1, l2) = fg.a;
%                     varN2(l, l1, l2) = fg.c;
%                     fg = fit(x.', real(s3(:, l, l1, l2)), renyi, 'StartPoint', [1 L * (v - u) / 2 0.1]);
%                     sn03(l, l1, l2) = fg.a;
%                     varN3(l, l1, l2) = fg.c;
%                     fg = fit(x.', real(s4(:, l, l1, l2)), renyi, 'StartPoint', [1 L * (v - u) / 2 0.1]);
%                     sn04(l, l1, l2) = fg.a;
%                     varN4(l, l1, l2) = fg.c;
%                     fg = fit(x.', real(s5(:, l, l1, l2)), renyi, 'StartPoint', [1 L * (v - u) / 2 0.1]);
%                     sn05(l, l1, l2) = fg.a;
%                     varN5(l, l1, l2) = fg.c;
                end
            end
        end
    end
    clear u;
    clear v;
    clear l;
    clear l1;
    clear l2;
    clear cicj;
    clear eigenVals;
    clear f;
    save(filename);
end
