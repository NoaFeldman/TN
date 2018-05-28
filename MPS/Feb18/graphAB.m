function graphAB(L)
    % recreate Fig 2 here https://arxiv.org/pdf/1711.09418.pdf
    n = 1; K = 1;
    x = (L/2 - 5):(L/2 + 5);
    avgNA = L / 2;
%     sigma2 = log(L) + 1 + 0.577 + log(2);
    sigma2 = log(L);
    pTheo = sqrt(pi * n / (2 * K * sigma2)) .* exp(-1 * n * pi^2 .* (x - avgNA).^2 ./ (2 * K * sigma2));
    plot(x, pTheo, 'k');
    hold on
    legendInfo{1} = ['$P(N_A)$'];
    STheo = sqrt(pi * sigma2 / (2 * K)) / 3 .* exp(-1 * pi^2 .* (x - avgNA).^2 ./ (2 * K * sigma2));
    plot(x, STheo, 'r');
    legendInfo{2} = ['$S(N_A)$'];
    
    C = zeros(L, L);
    for i = 1 : L
        for j = 1 : L
            if (i == j) 
                C(i, j) = 0.5; %2;
            else
                C(i, j) = sin(pi * (i-j) / 2) / (pi * (i - j));
                C(j, i) = C(i, j);
            end
        end
    end
    [~, f] = eig(C);
    e = zeros(1, L);
    for i = 1 : L
        e(i) = log((1 - f(i, i))/f(i, i));
    end
    p = getSNA(1, e, x);
    scatter(x(:), p(:), 'markerEdgeColor', 'k');
    
    S = getEE(e, x);
    scatter(x(:), S(:), 'markerEdgeColor', 'r');
        
    legend(legendInfo, 'Interpreter', 'latex');
    xlabel('$N_A$', 'Interpreter', 'latex');
%     savefig('FigAB');
    save(strcat('graphAB', int2str(L)), 'p', 'S');
    hold off
end

function S = getEE(e, x)
    dn = 1e-3;
    snTop = getSNA(1 + dn / 2, e, x);
    snBottom = getSNA(1 - dn / 2, e, x);
    S = -(snTop - snBottom) / dn;
end

function s = getSNA(n, e, x) 
    NSteps = 1e4;
    stepSize = 2 * pi / NSteps ;
    s = zeros(1, length(x));
    for alpha = -pi : stepSize : pi
        s1Alpha = getSAlpha(n, alpha, e);
        s = s + s1Alpha .* exp(complex(0, -alpha .* x)) / (2 * pi / stepSize);
    end
end

function s = getSAlpha(n, alpha, e)
    s = 1;
    f = 1 ./ (exp(e) + 1);
    for l = 1 : length(e)
        s = s * (exp(j * alpha) * f(l)^n + (1 - f(l))^n);
    end
end
