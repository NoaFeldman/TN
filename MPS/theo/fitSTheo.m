function fitSTheo(data, L, pointFunc, delta, cftRegion)
    % Fitting the numeric / exact results of the chargereolved entropies to
    % theory, based on:
    % s_n(\alpha) = (\gamma_1 L_eff)^{c pointFunc / 12 (n - n^-1)}
    % (\gamma_2 L_eff)^{pointFunc K / n (\alpha/2 \pi)^2}
    % s_n(N_A) = \int_-pi ^ pi s_n(\alpha) e^{-i N_A \alpha} d\alpha / 2\pi
    % without taking the integral boundaries to infinity.
    K = getK(delta);
    
    hold off;
    [as, ~, ~, ~, yfit] = fitnonlin(data.t(cftRegion), data.t(cftRegion), ...
        real(data.sFull(cftRegion)), 0.01.*data.t(cftRegion), 0.01.*real(data.sFull(cftRegion)), ...
        'getEntanglementEntropy', [1 1], [L 1]);
    plot(data.t(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(data.t(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(data.t, as(1), L, 1);
    hold off
    
    f = fittype('exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    sig2 = zeros(length(data.t), 1);
    for i = 1:length(data.t)
        [fg, gof] = fit(data.alphas.', real(data.s1Alpha(:, i)), f, 'StartPoint', 0.5);
        plot(fg, data.alphas, real(data.s1Alpha(:, i))); pause(0.01);
        sig2(i) = fg.c;
    end
    
    [ap, ~, ~, ~, yfit] = fitnonlin(data.t(cftRegion), data.t(cftRegion), ...
        sig2(cftRegion).', 0.01.*data.t(cftRegion), ...
        0.01.*sig2(cftRegion).', 'getSigmaAlpha', [1e-4 1], [L 1 K pointFunc]);
    plot(data.t(cftRegion), sig2(cftRegion));
    hold on
    plot(data.t(cftRegion), yfit);
    L2 = exp(pi^2./getSigmaAlpha(data.t, ap, [L 1 K pointFunc]));
    hold off
    
    plot(data.t, data.s(6, :), 'color', 'c');
    hold on
    plot(data.t, data.s(5, :), 'color', 'c');
    plot(data.t, data.s(4, :), 'color', 'c');
    max = length(data.t);
    plot(data.t(2:max), stheo(L1(2:max), L2(2:max), 0), 'color', 'b'); 
    plot(data.t(2:max), stheo(L1(2:max), L2(2:max), 1), 'color', 'b'); 
    plot(data.t(2:max), stheo(L1(2:max), L2(2:max), 2), 'color', 'b');   
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    title('Charged resolved entanglement entropy');
    h = findobj(gca);
    legend([h(5) h(2)], {'Exact', 'CFT'});
end

function K = getK(delta)
    K = pi ./ (2 * (pi - acos(delta)));
end