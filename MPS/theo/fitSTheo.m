function fitSTheo(sz, alphas, t, fitStartInd, fitEndInd, sigmaN, sFull, s, filename, figname, L)
    %
    % Find var(N_A) and fit to theoretic expectation to sigmaN.
%     sigmaN = zeros(fitEndInd, 1);
%     avgN = zeros(fitEndInd, 1);
%     R = zeros(fitEndInd, 1);
%     f = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
%     for i = fitStartInd:fitEndInd
%         [fg, gof] = fit(sz.', real(s1Charge(:, i)), f, 'StartPoint', [0.1 1]);
%         sigmaN(i) = fg.c;
%         avgN(i) = fg.b;
%         R(i) = gof.rsquare;
%     end
%     [ap, aerrp ,covp, chi2p, yfit] = fitnonlin(t(fitStartInd:fitEndInd), t(fitStartInd:fitEndInd), ...
%         sigmaN(fitStartInd:fitEndInd).', 0.01.*t(fitStartInd:fitEndInd), ...
%         0.01.*sigmaN(fitStartInd:fitEndInd).', 'getSigmaN', [15 1], L);
%     mysig = getSigmaN(t, ap, L);
%     epsP = ap(1);
%     constP = ap(2);
    
    % Fit flux-resolved Renyi entropies for the flux / charge dependent
    % parameters.
%     f = fittype('b^(-x^2) + w*b^(-(x + 2*pi)^2) + w*b^(-(x - 2*pi)^2)', 'independent', 'x', 'dependent', 'y');
%     XChargeScaled = zeros(1, length(t));
%     w = zeros(1, length(t));
%     for i = fitStartInd:fitEndInd
%         [fg, gof] = fit(alphas.', real(s1Alpha(:, i)), f, 'StartPoint', [1 0.85]);
% %         plot(fg, alphas.', real(s1Alpha(:, i))); pause(0.5);
%         XChargeScaled(i) = fg.b;
%         w(i) = fg.w;
%     end

    

    % Fit full (non-sectored) entanglement entropy to theoretical
    % expectation.
    [as, aerrs, covs, chi2s, yfit] = fitnonlin(t(fitStartInd:fitEndInd), t(fitStartInd:fitEndInd), ...
        real(sFull(fitStartInd:fitEndInd)), 0.01.*t(fitStartInd:fitEndInd), ...
        0.01.*real(sFull(fitStartInd:fitEndInd)), 'getEntanglementEntropy', [1 1], [L, 1]);
    epsS = as(1);
    c1 = as(2);
    X0Scaled = exp(6.*(sFull - 1/6.*log(getScaledVariable(t, epsS, L, 1)))).*getScaledVariable(t, epsS, L, 1);
    XChargeScaled = exp(2*pi^2.*sigmaN);
    % Use fit parameters for plotting theoretical expectation for sectored
    % entanglement entropy.
    hold off;
    plot(t, s(round(length(sz) / 2), :), 'color', 'c');
    hold on;
    plot(t, s(round(length(sz) / 2) + 1, :), 'color', 'c');
    plot(t, s(round(length(sz) / 2) + 2, :), 'color', 'c');
    plot(t(fitStartInd:fitEndInd), stheo(X0Scaled(fitStartInd:fitEndInd), XChargeScaled(fitStartInd:fitEndInd), 0, ...
        1, 1, 1), 'color', 'b');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    title('Charged resolved entanglement entropy');
    h = findobj(gca);
    legend([h(5) h(2)], {'Exact', 'CFT'});
end
    