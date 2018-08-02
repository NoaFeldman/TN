function analyzeSAlpha(filename)
    hold off
    % load('t0L10000_1point_full.mat');
    data = load(filename);
    L = 10000;
    % Plot sigma^2 extracted from s_{1-5}
    % From Eq.(7), we expect it to go as 1/n
    for n = 1:5
        plot(data.ratios, permute(n.*data.(strcat('varN', int2str(n))), [3 1 2]));
        legendinfo{n} = strcat('n = ', int2str(n));
        hold on;
    end
    hold off
       
    % Evaluating |S_1(alpha)|
    f = fittype('b^(-x^2)', 'independent', 'x', 'dependent', 'y');
    % We expact b \propto (L*sin(pi*l/L))^(k/n*(2pi)^2).
    b1s = zeros(1, length(data.ratios));
    fittingRange = 20:(length(data.alphas) - 20);
    for i = 1:length(data.ratios)
        [fg, gof] = fit(data.alphas(fittingRange).', real(data.s1Alpha(fittingRange, 1, 1, i)), f, ...
            'StartPoint', [(L*sin(pi*data.ratios(i)))^(1/(2*pi)^2)]);
        plot(fg, data.alphas.', abs(data.s1Alpha(:, 1, 1, i))); % pause(0.5);
        b1s(i) = fg.b;
    end
     
    % Note the corners it gets here near l = 0.5L!
    plot(data.ratios, b1s, 'color', 'r');
    hold on
    scatter(data.ratios, b1s, '+', 'markerEdgeColor', 'r');
    plot(data.ratios, (L*sin(pi.*data.ratios)).^(1/(2*pi)^2), 'color', 'b');
    hold off
    % Note we only expected this to be proportional (no pi, c1,
    % 2-factor...), isn't it wierd?
    plot(data.ratios, (L*sin(pi.*data.ratios)).^(1/(2*pi)^2)./b1s);
    hold off
    
    % TODO - go over S_{2,3}(\alpha)
end