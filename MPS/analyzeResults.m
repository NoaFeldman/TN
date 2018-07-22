function analyzeResults(var, func, sFull, sz, p, s, pointFunc)
    % Var is our x variable
    % func(var) is the scaled quantity
    % S_A = pointFunc/6*[log(func(var)) + c]
    hold off;
    
    scatter(log(func(var)), sFull);
    hl = lsline;
    % We expect sFullLogCoeff = pointFunc/6
    sFullLogCoeff = (hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2));
    sFullC1 = (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2));

    scaled = exp(-(6/pointFunc).*sFullC1).*func(var).^(-1);
    scaledEmpiric = exp(-(6/pointFunc).*(sFull - 1/(6/pointFunc).*log(func(var)))).*func(var).^(-1);
   
    sigmaNEmpiric = zeros(1, length(var));
    R = zeros(1, length(var));
    f = fittype('sqrt(1/(2*pi*c))*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    for i = 1:length(var)
        [fg, gof] = fit(sz.', real(p(:, i)), f, 'StartPoint', [1]); 
        sigmaNEmpiric(i) = fg.c;
        R(i) = gof.rsquare;
    end
    scatter(log(func(var)), sigmaNEmpiric);
    hl = lsline;
    % We expect sigmaLogCoeff = pointFunc/(2*pi^2)
    sigmaLogCoeff = (hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2));
    sigmaC1 = (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2));
    sigmaN = log(func(var))*sigmaLogCoeff + sigmaC1;
    
    plot(var, s(6, :), 'color', 'c');
    hold on
    plot(var, s(5, :), 'color', 'c');
    plot(var, s(4, :), 'color', 'c');
    plot(var, stheo(scaledEmpiric, sigmaNEmpiric, 0, pointFunc), 'color', 'b');
    plot(var, stheo(scaledEmpiric, sigmaNEmpiric, 1, pointFunc), 'color', 'b');
    plot(var, stheo(scaledEmpiric, sigmaNEmpiric, 2, pointFunc), 'color', 'b');
end