function fitSTheo(sz, t, fitStartInd, fitEndInd, p, sFull, s, filename, figname)
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
        0.01.*sigmaN(fitStartInd:fitEndInd).', 'getSigmaN', [25 1 1]);
    mysig = getSigmaN(t, ap);
    epsP = ap(1);
    constP = ap(3);
    % Fit full (non-sectored) entanglement entropy to theoretical
    % expectation.
    [as, aerrs, covs, chi2s, yfit] = fitnonlin(t(fitStartInd:fitEndInd), t(fitStartInd:fitEndInd), ...
        real(sFull(fitStartInd:fitEndInd)), 0.01.*t(fitStartInd:fitEndInd), ...
        0.01.*real(sFull(fitStartInd:fitEndInd)), 'getEntanglementEntropy', [1 1]);
    epsS = as(1);
    c1 = as(2);
    % Use fit parameters for plotting theoretical expectation for sectored
    % entanglement entropy.
    scatter(t, s(6, :), 'markerEdgeColor', 'c');
    hold on;
    scatter(t, s(7, :), 'markerEdgeColor', 'c');
    scatter(t, s(8, :), 'markerEdgeColor', 'c');
    plot(t, stheo(t, 0, epsS, c1, mysig), 'color', 'b');
    plot(t, stheo(t, 1, epsS, c1, mysig), 'color', 'b');
    plot(t, stheo(t, 2, epsS, c1, mysig), 'color', 'b');
    save(filename);
    savefig(figname);
end
    