function analyzeSAlpha(filename)
    hold off
    % load('t0L10000_1point_full.mat');
    data = load(filename);
    L = 10000;
    
    % Index for location of alpha = 0 in data.alphas.
    zeroIndex = 63;
    
    % Plot sigma^2 extracted from s_{1-5}
    % From Eq.(7), we expect it to go as 1/n
%     for n = 1:5
%         plot(data.ratios, permute(n.*data.(strcat('varN', int2str(n))), [3 1 2]));
%         legendinfo{n} = strcat('n = ', int2str(n));
%         hold on;
%     end
%     hold off
       
    % Evaluating S_1(alpha)
    f = fittype('b^(-x^2) + b^(-(x + 2*pi)^2) + b^(-(x - 2*pi)^2)', 'independent', 'x', 'dependent', 'y');
%     f = fittype('b^(-x^2)', 'independent', 'x', 'dependent', 'y');
    
    % We expact b \propto (L*sin(pi*l/L))^(k/n*(2pi)^2).
    scaled = zeros(1, length(data.ratios));
    fittingRange = 20:(length(data.alphas) - 20);
    for i = 1:length(data.ratios)
        [fg, gof] = fit(data.alphas(fittingRange).', real(data.s1Alpha(fittingRange, 1, 1, i)), f, ...
            'StartPoint', [(L*sin(pi*data.ratios(i)))^(1/(2*pi)^2)]);
%         plot(fg, data.alphas.', abs(data.s1Alpha(:, 1, 1, i))); pause(0.5);
        scaled(i) = fg.b;
    end
     
    % Compare to the equivalent of GS Eq(7) with the new form for f.
    n = 1;
    s0 = 1;
    for dNA = 0:2
        plot(data.ratios, permute(data.s1(6 - dNA, 1, 1, :), [4 1 2 3]), 'color', 'c');
        hold on
        plot(data.ratios, getNewSNA(s0, n, dNA, scaled), 'color', 'b')
    end
    hold off

%     plot(data.ratios, b1s, 'color', 'r');
%     hold on
%     scatter(data.ratios, b1s, '+', 'markerEdgeColor', 'r');
%     plot(data.ratios, (L*sin(pi.*data.ratios)).^(1/(2*pi)^2), 'color', 'b');
%     hold off
%     % Note we only expected this to be proportional (no pi, c1,
%     % 2-factor...), isn't it wierd?
%     plot(data.ratios, (L*sin(pi.*data.ratios)).^(1/(2*pi)^2)./b1s);
%     hold off
    
    % Fit results for larger n values
    for n = 2:5
        s0 = zeros(1, length(data.ratios));
        for i = 1:length(data.ratios)
            sAlpha = data.(strcat('s', int2str(n), 'Alpha'));
            s0(i) = sAlpha(zeroIndex, 1, 1, i)./...
                (scaled(i).^(-data.alphas(zeroIndex).^2/n) + ...
                    scaled(i).^(-(data.alphas(zeroIndex) + 2*pi).^2/n) + ...
                    scaled(i).^(-(data.alphas(zeroIndex) - 2*pi).^2/n));
            hold off
            scatter(data.alphas, real(sAlpha(:, 1, 1, i)), '.', 'markerEdgeColor', 'b');
            hold on
            plot(data.alphas, s0(i).*...
                (scaled(i).^(-data.alphas.^2/n) + ...
                    scaled(i).^(-(data.alphas + 2*pi).^2/n) + ...
                    scaled(i).^(-(data.alphas - 2*pi).^2/n)), 'color', 'r'); 
            title(strcat('n = ', int2str(n)));
%             pause(0.5)
        end
        hold off
        for dNA = 0:2
            s = data.(strcat('s', int2str(n)));
            plot(data.ratios, permute(s(6 - dNA, 1, 1, :), [4 1 2 3]), 'color', 'c');
            hold on
            plot(data.ratios, getNewSNA(s0, n, dNA, scaled), 'color', 'b')
        end
    hold off
    end
end

function sna = getNewSNA(s0, n, dNA, scaled) 
    sna = s0.*exp(-dNA^2./4.*n./log(scaled))*1./(4*sqrt(pi)).*sqrt(n./log(scaled)).* ...
            (erfz(3*pi.*sqrt(log(scaled)./n) + 1i*dNA./2 .*sqrt(n./log(scaled))) - ...
                erfz(-3*pi.*sqrt(log(scaled)./n) + 1i*dNA/2 .*sqrt(n./log(scaled))));
end