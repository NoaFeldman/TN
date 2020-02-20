function fitSTheo(data, varName, L, pointFunc, K, u, model, cftRegion, figName)
    % Fitting the numeric / exact results of the chargereolved entropies to
    % theory, based on:
    % s_n(\alpha) = (\gamma_1 L_eff)^{c pointFunc / 12 (n - n^-1)}
    % (\gamma_2 L_eff)^{pointFunc K / n (\alpha/2 \pi)^2}
    % s_n(N_A) = \int_-pi ^ pi s_n(\alpha) e^{-i N_A \alpha} d\alpha / 2\pi
    % without taking the integral boundaries to infinity.
    
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
        'getEntanglementEntropy', [1e-4 0.2], [L model pointFunc u], []);
    fig = figure('visible','off');
    plot(var(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(var(cftRegion), yfit);
    L1 = exp((6/pointFunc).*as(2)).*getScaledVariable(var, as(1), L, model, u);
    saveas(fig, ['sFullFit' num2str(K)],'png');
    hold off
    
    f = fittype('w0* exp(-(x)^2 / (2 * s)) +  w1 * (exp(-(x - 2*pi)^2 / (2 * s)) + exp(-(x + 2*pi)^2 / (2 * s)))', 'independent', 'x', 'dependent', 'y');
    V = zeros(1, length(var));
    w0s = zeros(1, length(var));
    w1s = zeros(1, length(var));
    alphaRegion = 1:length(data.alphas) - 1 + 1;
    for i = 1:length(var)
        [fg, gof] = fit(data.alphas(alphaRegion).', real(data.s1Alpha(alphaRegion, i)), ...
            f, 'StartPoint', [0.5 1 1]);
        if (mod(i, 1000) == 0)
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
        w0s(i) = fg.w0;
        w1s(i) = fg.w1;
    end
  
    w0 = mean(w0s(cftRegion));
    w1 = mean(w1s(cftRegion));
    
    [ap, ~, ~, ~, yfit] = fitnonlin(var(cftRegion), var(cftRegion), ...
        V(cftRegion).', 0.01.*var(cftRegion), ...
        0.01.*V(cftRegion).', 'getSigmaAlpha', [1e-2 2], [L model K pointFunc u], []);
    plot(var(cftRegion), V(cftRegion));
    hold on
    plot(var(cftRegion), yfit);
    L2 = exp(pi^2./(K .* getSigmaAlpha(var, ap, [L model K pointFunc u])));
    saveas(fig, ['sigmaAlphaFit' num2str(K)],'png');    
    hold off
    
    scatter(var, data.sFull, '.', 'markerEdgeColor', [0 0 0.3]);
    hold on
    scatter(var, data.s(6, :), '.', 'markerEdgeColor', [0 0 0.3]);
    scatter(var, data.s(5, :), '.', 'markerEdgeColor', [0 0 0.3]);
    scatter(var, data.s(4, :), '.', 'markerEdgeColor', [0 0 0.3]);
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 0, K, 1, 0), '--', 'color', [0.9 0 0.4]);
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 1, K, 1, 0), '--', 'color', [0.9 0 0.4]);
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 2, K, 1, 0), '--', 'color', [0.9 0 0.4]);
    plot(var(3:end), getEntanglementEntropy(var(3:end), as, [L model pointFunc u]), 'color', [0 0.9 0.4]);
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 0, K, w0, w1), 'color', [0 0.9 0.4]);
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 1, K, w0, w1), 'color', [0 0.9 0.4]); 
    plot(var(3:end), stheo2(L1(3:end), L2(3:end), 2, K, w0, w1), 'color', [0 0.9 0.4]);
    disp('CFT prediction too high at the beginning? Heere are the predictions so I know when to start:');
    disp(stheo2(L1(2:11), L2(2:11), 1, K, w0, w1));
    xlabel('$Jt$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$S_A(N_A)$', 'Interpreter', 'latex', 'FontSize', 20);
    h = findobj(gca);
    legend([h(8) h(2) h(5)], {'Exact', 'CFT($w_1 \ne 0$)', 'CFT($w_1 = 0$)'}, 'Interpreter', 'latex', 'FontSize', 20);
    xlim([min(data.t) max(data.t)]);
    annotation('textbox',[.2 .2 .2 .6],'String', 'Total','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .4],'String', '$|\Delta N_A| = 0$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .2],'String', '$|\Delta N_A| = 1$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    annotation('textbox',[.2 .2 .2 .0],'String', '$|\Delta N_A| = 2$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
    box on
    savefig(figName);
    saveas(fig, strcat(figName, '.png'));
end
