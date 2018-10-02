function fitSTheo(data, L, pointFunc, delta, cftRegion, figName)
    % Fitting the numeric / exact results of the chargereolved entropies to
    % theory, based on:
    % s_n(\alpha) = (\gamma_1 L_eff)^{c pointFunc / 12 (n - n^-1)}
    % (\gamma_2 L_eff)^{pointFunc K / n (\alpha/2 \pi)^2}
    % s_n(N_A) = \int_-pi ^ pi s_n(\alpha) e^{-i N_A \alpha} d\alpha / 2\pi
    % without taking the integral boundaries to infinity.
    K = getK(delta);
    
    counter = 1;
    for i = 2:length(data.t) - 1
        if (data.s(6, i) < data.s(6, i - 1) && data.s(6, i) < data.s(6, i + 1))
            mins(counter) = i;
            counter = counter + 1;
        end
    end
    cftRegion = mins;
    plot(data.t, data.sFull)
    hold on
    scatter(data.t(cftRegion), data.sFull(cftRegion))
    
    hold off;
    [as, ~, ~, ~, yfit] = fitnonlin(data.t(cftRegion), data.t(cftRegion), ...
        real(data.sFull(cftRegion)), 0.01.*data.t(cftRegion), 0.01.*real(data.sFull(cftRegion)), ...
        'getEntanglementEntropy', [1 1], [L 1]);
    plot(data.t(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(data.t(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(data.t, as(1), L, 1);
    hold off
    
    f = fittype('exp(-x^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    var = zeros(1, length(data.t));
    f2 = fittype('b* exp(-(x)^2/(2*s)) + d * exp(-(x - 2*pi)^2/(2*s)) + d * exp(-(x + 2*pi)^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    var2 = zeros(1, length(data.t));
    w1 = zeros(1, length(data.t));
    w2 = zeros(1, length(data.t));
    alphaRegion = 50:580;
    for i = 1:length(data.t)
        [fg, gof] = fit(data.alphas(alphaRegion).', real(data.s1Alpha(alphaRegion, i)), ...
            f, 'StartPoint', [0.5]);
        [fg2, gof] = fit(data.alphas.', real(data.s1Alpha(:, i)), ...
            f2, 'StartPoint', [1 -1 0.5]);
        if (mod(i, 10) == 0)
%             plot(fg, data.alphas, real(data.s1Alpha(:, i))); 
%             hold on
%             plot(fg2, data.alphas, real(data.s1Alpha(:, i))); pause(0.01);
%             hold off
        end
        var(i) = fg.s;
        var2(i) = fg2.s;
        w1(i) = fg2.b;
        w2(i) = fg2.d;
    end
  
    w1 = mean(w1);
    w2 = mean(w2);
    
    [ap, ~, ~, ~, yfit] = fitnonlin(data.t(cftRegion), data.t(cftRegion), ...
        var2(cftRegion), 0.01.*data.t(cftRegion), ...
        0.01.*var2(cftRegion), 'getSigmaAlpha', [1e-4 1], [L 1 K pointFunc]);
    plot(data.t(cftRegion), var2(cftRegion));
    hold on
    plot(data.t(cftRegion), yfit);
    L2 = exp(pi^2./(K .* getSigmaAlpha(data.t, ap, [L 1 K pointFunc])));
%     L2 = exp(pi^2./(K.*var2));
    hold off
    
    plot(data.t, data.s(6, :), 'color', 'c');
    hold on
    plot(data.t, data.s(5, :), 'color', 'c');
    plot(data.t, data.s(4, :), 'color', 'c');
    max = length(data.t);
%     plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 0, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 1, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 2, K, w1(2:max), w2(2:max)), 'color', 'b');     
    plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 0, K, w1, w2), 'color', 'b'); 
    plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 1, K, w1, w2), 'color', 'b'); 
    plot(data.t(2:max), stheo2(L1(2:max), L2(2:max), 2, K, w1, w2), 'color', 'b');     
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    title(strcat('Charged resolved entanglement entropy, $\Delta = ', num2str(delta), '$'), ...
        'Interpreter', 'latex');
    h = findobj(gca);
    legend([h(5) h(2)], {'DMRG', 'CFT'});
    saveas(gcf, figName);
end
