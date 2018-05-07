function fitSTheo(sz, t, fitStartInd, fitEndInd, p, sFull, s, filename, figname, L)
    %
    % Find var(N_A) and fit to theoretic expectation to sigmaN.
    sigmaN = zeros(fitEndInd, 1);
    avgN = zeros(fitEndInd, 1);
    R = zeros(fitEndInd, 1);
    f = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    for i = fitStartInd:fitEndInd
        [fg, gof] = fit(sz.', real(p(:, i)), f, 'StartPoint', [0.1 1]);
        sigmaN(i) = fg.c;
        avgN(i) = fg.b;
        R(i) = gof.rsquare;
    end
    [ap, aerrp ,covp, chi2p, yfit] = fitnonlin(t(fitStartInd:fitEndInd), t(fitStartInd:fitEndInd), ...
        sigmaN(fitStartInd:fitEndInd).', 0.01.*t(fitStartInd:fitEndInd), ...
        0.01.*sigmaN(fitStartInd:fitEndInd).', 'getSigmaN', [25 1], L);
    mysig = getSigmaN(t, ap, L);
    epsP = ap(1);
    constP = ap(2);
    % Fit full (non-sectored) entanglement entropy to theoretical
    % expectation.
    [as, aerrs, covs, chi2s, yfit] = fitnonlin(t(fitStartInd:fitEndInd), t(fitStartInd:fitEndInd), ...
        real(sFull(fitStartInd:fitEndInd)), 0.01.*t(fitStartInd:fitEndInd), ...
        0.01.*real(sFull(fitStartInd:fitEndInd)), 'getEntanglementEntropy', [1 1], L);
    epsS = as(1);
    c1 = as(2);
    % Use fit parameters for plotting theoretical expectation for sectored
    % entanglement entropy.
%     hold off;
    plot(t, s(round(length(sz) / 2), :), 'color', 'c');
    hold on;
    plot(t, s(round(length(sz) / 2) + 1, :), 'color', 'c');
    plot(t, s(round(length(sz) / 2) + 2, :), 'color', 'c');
    plot(t, stheo(t, 0, epsS, c1, mysig, L), 'color', 'm');
    plot(t, stheo(t, 1, epsS, c1, mysig, L), 'color', 'm');
    plot(t, stheo(t, 2, epsS, c1, mysig, L), 'color', 'm');
%     plot(t, stheo(t, 0, epsS, c1, mysig, L), 'color', 'k');
%     plot(t, stheo(t, 1, epsS, c1, mysig, L), 'color', 'k');
%     plot(t, stheo(t, 2, epsS, c1, mysig, L), 'color', 'k');
    save(filename);
    savefig(figname);
end
    