function fitSTheo(data, varName, L, pointFunc, delta, model, cftRegion, figName)
    % Fitting the numeric / exact results of the chargereolved entropies to
    % theory, based on:
    % s_n(\alpha) = (\gamma_1 L_eff)^{c pointFunc / 12 (n - n^-1)}
    % (\gamma_2 L_eff)^{pointFunc K / n (\alpha/2 \pi)^2}
    % s_n(N_A) = \int_-pi ^ pi s_n(\alpha) e^{-i N_A \alpha} d\alpha / 2\pi
    % without taking the integral boundaries to infinity.
    K = getK(delta);
    
    var = data.(varName);
    
%     counter = 1;
%     for i = 2:length(var) - 1
%         if (data.s(6, i) < data.s(6, i - 1) && data.s(6, i) < data.s(6, i + 1))
%             mins(counter) = i;
%             counter = counter + 1;
%         end
%     end
%     cftRegion = mins;
%     plot(var, data.sFull)
%     hold on
%     scatter(var(cftRegion), data.sFull(cftRegion))
    
    hold off;
    [as, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        real(data.sFull(cftRegion)).', 0.01.*var(cftRegion), 0.01.*real(data.sFull(cftRegion)).', ...
        'getEntanglementEntropy', [1 1], [L model pointFunc]);
    plot(var(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(var(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(var, as(1), L, model);
    hold off
    
%     f = fittype('exp(-x^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    V = zeros(1, length(var));
    f2 = fittype('exp(-(x)^2/(2*s)) +  2 * 1 * exp(-((x - 2*pi)^2 + (x + 2*pi)^2)/(2*2*s))', 'independent', 'x', 'dependent', 'y');
%     f2 = fittype('exp(-(x)^2/(2*s)) + 1e-4 * exp(-(x - 2*pi)^2/(2*s)) + 1e-4 * exp(-(x + 2*pi)^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    V2 = zeros(1, length(var));
    w1s = zeros(1, length(var));
    w2s = zeros(1, length(var));
    alphaRegion = 1:125;
    for i = 1:length(var)
%         [fg, gof] = fit(data.alphas(alphaRegion).', real(data.s1Alpha(alphaRegion, i)), ...
%             f, 'StartPoint', [0.5]);
        [fg2, gof] = fit(data.alphas.', real(data.s1Alpha(:, i)), ...
            f2, 'StartPoint', [0.5]);
        if (mod(i, 1) == 0)
%             plot(fg, data.alphas, real(data.s1Alpha(:, i))); 
%             hold on
            plot(fg2, data.alphas, real(data.s1Alpha(:, i))); pause(0.1);
            hold off
        end
%         V(i) = fg.s;
        V2(i) = fg2.s;
%         w1s(i) = fg2.w1;
%         w2s(i) = fg2.w2;
    end
  
    w1 = 1; %mean(w1s(cftRegion));
    w2 = 1.5; %mean(w2s(cftRegion));
    
    [ap, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        V2(cftRegion).', 0.01.*var(cftRegion), ...
        0.01.*V2(cftRegion).', 'getSigmaAlpha', [1e-4 1], [L model K pointFunc]);
    plot(var(cftRegion), V2(cftRegion));
    hold on
    plot(var(cftRegion), yfit);
    L2 = exp(pi^2./(K .* getSigmaAlpha(var, ap, [L model K pointFunc])));
%     L2 = exp(pi^2./(K.*var2));
    hold off
    
    plot(var, data.s(6, :), 'color', 'c');
    hold on
    plot(var, data.s(5, :), 'color', 'c');
    plot(var, data.s(4, :), 'color', 'c');
    max = length(var);
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 0, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 1, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 2, K, w1(2:max), w2(2:max)), 'color', 'b');     
    plot(var(2:max), stheo2(L1(2:max), L2(2:max), 0, K, w1, w2), 'color', 'b'); 
    plot(var(2:max), stheo2(L1(2:max), L2(2:max), 1, K, w1, w2), 'color', 'b'); 
    plot(var(2:max), stheo2(L1(2:max), L2(2:max), 2, K, w1, w2), 'color', 'b');     
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    title(strcat('Charged resolved entanglement entropy, $\Delta = ', num2str(delta), '$'), ...
        'Interpreter', 'latex');
    h = findobj(gca);
    legend([h(5) h(2)], {'DMRG', 'CFT'});
    saveas(gcf, figName);
end
