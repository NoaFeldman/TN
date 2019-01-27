function fitSTheo(data, varName, L, pointFunc, delta, model, cftRegion, figName)
    % Fitting the numeric / exact results of the chargereolved entropies to
    % theory, based on:
    % s_n(\alpha) = (\gamma_1 L_eff)^{c pointFunc / 12 (n - n^-1)}
    % (\gamma_2 L_eff)^{pointFunc K / n (\alpha/2 \pi)^2}
    % s_n(N_A) = \int_-pi ^ pi s_n(\alpha) e^{-i N_A \alpha} d\alpha / 2\pi
    % without taking the integral boundaries to infinity.
    K = getK(delta);
    
    var = data.(varName);
        
    hold off;
    epsilon = 1e-4;
%     logLeff = log(getScaledVariable(var, epsilon, L, model));
%     scatter(logLeff(cftRegion), real(data.sFull(cftRegion)));
%     hold on
%     hl = lsline;
%     B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), ...
%         (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
%     L1 = exp((1/(pointFunc*B(1))).*B(2)).*getScaledVariable(var, epsilon, L, model);
    [as, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        real(data.sFull(cftRegion)), 0.01.*var(cftRegion), 0.01.*real(data.sFull(cftRegion)), ...
        'getEntanglementEntropy', [1e-4 0.2], [L model pointFunc], []);
    plot(var(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(var(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(var, as(1), L, model);
    hold off
    
    f = fittype('w1* exp(-(x)^2 / (2 * s)) +  w2 * (exp(-(x - 2*pi)^2 / (2 * s)) + exp(-(x + 2*pi)^2 / (2 * s)))', 'independent', 'x', 'dependent', 'y');
%     f = fittype('exp(-(x)^2 / (2 * s))', 'independent', 'x', 'dependent', 'y');
    V = zeros(1, length(var));
    w1s = zeros(1, length(var));
    w2s = zeros(1, length(var));
    alphaRegion = 1:length(data.alphas) - 1 + 1;
    for i = 1:length(var)
        [fg, gof] = fit(data.alphas(alphaRegion).', real(data.s1Alpha(alphaRegion, i)), ...
            f, 'StartPoint', [0.5 1 1]);
%         fFixed = fittype(strcat('w1* exp(-(x)^2 / (2 * s)) +  ', num2str(fg.w2), ' * (exp(-(x - 2*pi)^2 / (2 * s)) + exp(-(x + 2*pi)^2 / (2 * s)))'), ...
%             'independent', 'x', 'dependent', 'y');
%         [fg2, gof] = fit(data.alphas(alphaRegion).', real(data.s2Alpha(alphaRegion, i)./data.s2Alpha(315, i)), ...
%             fFixed, 'StartPoint', [0.5 1]);
%         [fg3, gof] = fit(data.alphas(alphaRegion).', real(data.s3Alpha(alphaRegion, i)./data.s3Alpha(315, i)), ...
%             fFixed, 'StartPoint', [0.5 1]);
%         [fg4, gof] = fit(data.alphas(alphaRegion).', real(data.s4Alpha(alphaRegion, i)./data.s4Alpha(315, i)), ...
%             fFixed, 'StartPoint', [0.5 1]);
%         [fg5, gof] = fit(data.alphas(alphaRegion).', real(data.s5Alpha(alphaRegion, i)./data.s5Alpha(315, i)), ...
%             fFixed, 'StartPoint', [0.5 1]);
        if (mod(i, 300) == 0)
%             plot(fg, data.alphas, real(data.s1Alpha(:, i))); pause(0.1);
            scatter(data.alphas, real(data.s1Alpha(:, i)), '.', 'markeredgecolor', [0 0 0.3]);
            hold on
            f2 = fittype('exp(-(x)^2 / (2 * s))', 'independent', 'x', 'dependent', 'y');
            [fg2, gof] = fit(data.alphas.', real(data.s1Alpha(:, i)), ...
                f2, 'StartPoint', [0.5]);
            plot(data.alphas, fg2(data.alphas), 'color', [0.7 0 0.3])
            plot(data.alphas, fg(data.alphas), 'color', [0 0.9 0.4])
            xlabel('$\alpha$', 'Interpreter', 'latex');
            ylabel('$s_1(\alpha)$', 'Interpreter', 'latex');
            title(strcat('$L/a = ', num2str(L), ', t = ', num2str(var(i)), '$'), 'Interpreter', 'latex');
            legend({'Exact\ndiagonalization', '$w_\pm = 0$', '$w_\pm \ne 0$'}, 'interpreter', 'latex')
            rectangle('Position',[-pi-0.2 0 0.6 0.1],'Curvature', [1, 1], 'EdgeColor', [0 0.8 0.8])
            hold off
        end
        V(i) = fg.s;
%         V2(i) = fg2.s;
%         V3(i) = fg3.s;
%         V4(i) = fg4.s;
%         V5(i) = fg5.s;
        w1s(i) = fg.w1;
        ws(i) = fg.w2;
%         w2s(i) = fg2.w2;
%         w3s(i) = fg3.w2;
%         w4s(i) = fg4.w2;
%         w5s(i) = fg5.w2;
    end
  
    w1 = mean(w1s(cftRegion));
    w2 = mean(ws(cftRegion));
    
    [ap, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        V(cftRegion).', 0.01.*var(cftRegion), ...
        0.01.*V(cftRegion).', 'getSigmaAlpha', [1e-4 2], [L model K pointFunc], []);
    plot(var(cftRegion), V(cftRegion));
    hold on
    plot(var(cftRegion), yfit);
    L2 = exp(pi^2./(K .* getSigmaAlpha(var, ap, [L model K pointFunc])));
    hold off
    w2
    ap(2)
    
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
    xlabel('$t[\hbar/J]$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$S(N_A)$', 'Interpreter', 'latex', 'FontSize', 16);
    title(strcat('Charged resolved entanglement entropy, $\Delta = ', num2str(delta), '$'), ...%, $w_\pm = 0$'), ...
        'Interpreter', 'latex', 'FontSize', 16);
    h = findobj(gca);
    legend([h(5) h(2)], {'Exact diagonalization', 'CFT'}, 'Interpreter', 'latex', 'FontSize', 14);
%     legend([h(5) h(2)], {'DMRG', 'CFT'}, 'Interpreter', 'latex', 'FontSize', 14);
    xlim([min(data.t) max(data.t)]);
    annotation('textbox',[.2 .2 .2 .6],'String', '$\Delta N_A = 0$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .4],'String', '$\Delta N_A = 1$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .2],'String', '$\Delta N_A = 2$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    box on
    savefig(figName);
    saveas(gcf, strcat(figName, '.png'));
end
