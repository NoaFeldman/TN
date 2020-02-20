function analyzeSAlpha(filename, L, cftRegion, varFieldName, model, pointFunc)
    hold off
    data = load(filename);
    
    zeroIndex = find(data.alphas == 0);

    var = data.(varFieldName);
    [~, varIndex] = max(data.sFull(:));
    [a, ~, ~, chi2, ~] = fitnonlin(data.alphas, ...
        data.alphas, real(data.s1Alpha(:, varIndex)).', ...
        0.01.*data.alphas, 0.01.*real(data.s1Alpha(:, varIndex)).', ...
        'sAlphaTheo', [-1e-4 2 0.95 5e9], [1 L var(varIndex) model]);
    plot(data.alphas, real(data.s1Alpha(:, varIndex)).');
    hold on
    plot(data.alphas, sAlphaTheo(data.alphas, a, [1 L var(varIndex) model]));
    
    for i = 10:10:length(var)
        hold off
        plot(data.alphas, real(data.s1Alpha(:, i)));
        hold on
        plot(data.alphas, sAlphaTheo(data.alphas, a, [1 L var(i) model]));
        pause(0.5);
    end
    
    hold off
    
    epsilon = a(1);
    multiplyingConst = a(2);
    w = a(3);
    XChargeScaled = multiplyingConst .* getScaledVariable(var, epsilon, L, model);
       
    
    [as, ~, ~, chi2, ~] = fitnonlin(var(cftRegion), ...
        var(cftRegion), real(data.sFull(cftRegion)), 0.01.*var(cftRegion), ...
        0.01.*real(data.sFull(cftRegion)), 'getEntanglementEntropy', [1 1], [L, model]);
    plot(var(cftRegion), real(data.sFull(cftRegion)));
    hold on
    plot(var(cftRegion), getEntanglementEntropy(var(cftRegion), as, [L, model]));
    
    hold off
    
    sFullEpsilon = as(1);
    sFullC1 = as(2);
    
    X0Scaled = exp((6/pointFunc).*sFullC1).*getScaledVariable(var, sFullEpsilon, L, model);
    
    plot(var, data.s(6, :), 'color', 'c');
    hold on
    plot(var, data.s(5, :), 'color', 'c');
    plot(var, data.s(4, :), 'color', 'c');
    k = 1;
    % Still need to play a little with multiplyingConst here - s1Alpha is
    % very insensitive to it but stheo is.
    plot(var(cftRegion), stheo(X0Scaled(cftRegion), XChargeScaled(cftRegion), 0, pointFunc, w, k), 'color', 'b');
    plot(var(cftRegion), stheo(X0Scaled(cftRegion), XChargeScaled(cftRegion), 1, pointFunc, w, k), 'color', 'b');
    plot(var(cftRegion), stheo(X0Scaled(cftRegion), XChargeScaled(cftRegion), 2, pointFunc, w, k), 'color', 'b');
end
