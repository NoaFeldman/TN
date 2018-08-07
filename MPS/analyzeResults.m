function analyzeResults(var, func, sFull, sz, p, s, sigmaNEmpiric, pointFunc)
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
   
%     sigmaNEmpiric = zeros(1, length(var));
%     R = zeros(1, length(var));
% %     f = fittype('sqrt(1/(2*pi*c))*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
%     f = fittype('b*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
%     for i = 1:length(var)
%         [fg, gof] = fit(sz.', real(p(:, i)), f, 'StartPoint', [1 0.1]); 
%         sigmaNEmpiric(i) = fg.c;
%         R(i) = gof.rsquare;
%     end
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
    hold off
    plot(var, (s(6, :) - stheo(scaledEmpiric, sigmaNEmpiric, 0, pointFunc)) ./ s(6, :), 'color', [0, 0.1, 1]);
    hold on;
    plot(var, (s(5, :) - stheo(scaledEmpiric, sigmaNEmpiric, 1, pointFunc)) ./ s(5, :), 'color', [0, 0.4470, 0.7410]);
    plot(var, (s(4, :) - stheo(scaledEmpiric, sigmaNEmpiric, 2, pointFunc)) ./ s(4, :), 'color', [0.3010, 0.7450, 0.9330]);
    plot(var, (sFull - sFullLogCoeff.*log(func(var)) - sFullC1) ./ sFull, 'color', 'r');
    legend({'$\Delta N_A = 0$', '$\Delta N_A = 1$', '$\Delta N_A = 2$', 'sFull'}, 'Interpreter', 'latex');
end