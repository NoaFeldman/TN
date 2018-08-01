function analyzeSAlpha()
    hold off
    load('t0L10000_1point_full.mat');
    % Plot sigma^2 extracted from s_{1-5}
    % From Eq.(7), we expect it to go as 1/n
    for n = 1:5
        plot(f10.ratios, n.*f10.(strcat('sigma', int2str(n))));
        legendinfo{n} = strcat('n = ', int2str(n));
        hold on;
    end
    hold off
       
    % Evaluating |S_1(alpha)|
    f = fittype('b^(-c*x^2)', 'independent', 'x', 'dependent', 'y');
    % We expact b \propto L*sin(pi*l/L), c = K/(2pi)^2.
    bs = zeros(1, length(f10.ratios));
    cs = zeros(1, length(f10.ratios));
    for i = 1:length(f10.ratios)
        [fg, gof] = fit(f10.alphas.', abs(f10.sAlpha(:, 1, 1, i)), f, ...
            'StartPoint', [10000*sin(pi*f10.ratios(i)) 1/(2*pi)^2]);
        plot(fg, f10.alphas.', abs(f10.sAlpha(:, 1, 1, i))); pause(0.5);
        bs(i) = fg.b;
        cs(i) = fg.c;
    end
    plot(f10.ratios, bs);
    hold on
    plot(f10.ratios, 10000*sin(pi.*f10.ratios));
    hold off
    plot(f10.ratios, cs);
end