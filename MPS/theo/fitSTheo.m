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
        real(data.sFull(cftRegion)), 0.01.*var(cftRegion), 0.01.*real(data.sFull(cftRegion)), ...
        'getEntanglementEntropy', [1 1], [L model pointFunc], []);
    plot(var(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(var(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(var, as(1), L, model);
    hold off
    
%     f = fittype('exp(-x^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    V = zeros(1, length(var));
    f2 = fittype('exp(-(x)^2/(2*s)) +  w2 * (exp(-(x - 2*pi)^2/(2*s)) + exp(-(x + 2*pi)^2/(2*s)))', 'independent', 'x', 'dependent', 'y');
%     f2 = fittype('exp(-(x)^2/(2*s)) + 1e-4 * exp(-(x - 2*pi)^2/(2*s)) + 1e-4 * exp(-(x + 2*pi)^2/(2*s))', 'independent', 'x', 'dependent', 'y');
    V2 = zeros(1, length(var));
    w1s = zeros(1, length(var));
    w2s = zeros(1, length(var));
    alphaRegion = 1:length(data.alphas);
    for i = 1:length(var)
%         [fg, gof] = fit(data.alphas(alphaRegion).', real(data.s1Alpha(alphaRegion, i)), ...
%             f, 'StartPoint', [0.5]);
        [fg2, gof] = fit(data.alphas.', real(data.s1Alpha(:, i)), ...
            f2, 'StartPoint', [0.5 0.8]);
        if (mod(i, 1) == 0)
%             plot(fg, data.alphas, real(data.s1Alpha(:, i))); 
%             hold on
            plot(fg2, data.alphas, real(data.s1Alpha(:, i))); pause(0.01);
            hold off
        end
%         V(i) = fg.s;
        V2(i) = fg2.s;
%         w1s(i) = fg2.w1;
        w2s(i) = fg2.w2;
    end
  
    w1 = 1; %mean(w1s(cftRegion));
    w2 = 1.5; %mean(w2s(cftRegion));
    
    [ap, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        V2(cftRegion).', 0.01.*var(cftRegion), ...
        0.01.*V2(cftRegion).', 'getSigmaAlpha', [1e-4 1], [L model K pointFunc], []);
    plot(var(cftRegion), V2(cftRegion));
    hold on
    plot(var(cftRegion), yfit);
    L2 = exp(pi^2./(K .* getSigmaAlpha(var, ap, [L model K pointFunc])));
    hold off
    
    scatter(var, data.s(6, :), '.', 'markerEdgeColor', [0 0 0.3]);
    hold on
    scatter(var, data.s(5, :), '.', 'markerEdgeColor', [0 0 0.3]);
    scatter(var, data.s(4, :), '.', 'markerEdgeColor', [0 0 0.3]);
    maxInd = length(var);
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 0, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 1, K, w1(2:max), w2(2:max)), 'color', 'b'); 
%     plot(var(2:max), stheo2(L1(2:max), L2(2:max), 2, K, w1(2:max), w2(2:max)), 'color', 'b');     
    plot(var(2:maxInd), stheo2(L1(2:maxInd), L2(2:maxInd), 0, K, w1, w2), 'color', [0 0.9 0.4]); 
    plot(var(2:maxInd), stheo2(L1(2:maxInd), L2(2:maxInd), 1, K, w1, w2), 'color', [0 0.9 0.4]); 
    plot(var(2:maxInd), stheo2(L1(2:maxInd), L2(2:maxInd), 2, K, w1, w2), 'color', [0 0.9 0.4]);     
    xlabel('$t[\hbar/J]$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    title(strcat('Charged resolved entanglement entropy, $\Delta = ', num2str(delta), '$'), ...
        'Interpreter', 'latex');
    h = findobj(gca);
    legend([h(5) h(2)], {'DMRG', 'CFT'}, 'Interpreter', 'latex', 'FontSize', 14);
    xlim([min(data.t) max(data.t)]);
    annotation('textbox',[.2 .2 .2 .6],'String', '$\Delta N_A = 0$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .4],'String', '$\Delta N_A = 1$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .2],'String', '$\Delta N_A = 2$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    box on
    saveas(gcf, figName);
end
