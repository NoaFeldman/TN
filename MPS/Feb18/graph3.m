function graph3()
    L = [100 200 400 600 800 1000 1500];
    filenames = ["spectrumN2000D1_1024.mat", "spectrumN2000D1_1024.mat", "spectrumN2000D_1024.mat", "spectrumN1500D1_1024.mat"];
    Delta = 1;
    for l = 1 : length(L)
        s = load(strcat("spectrumN", int2str(L(l)), "D", num2str(Delta), "_1024.mat"));
        marr = s.spectrum.keys();
        sarr = zeros(0);
        for m = 1 : length(marr)
            sarr = cat(1, sarr, s.spectrum(marr{m}));
        end
        sarr = sort(sarr);
        % sarr is sorted smallest to largest
        for i = 1:length(sarr)
            if (sarr(i) > 1e-9)
                sarr = sarr(i : length(sarr));
                break;
            end
        end
        lambdaMax = sarr(length(sarr));
        b = - log(lambdaMax);
        x = zeros(length(sarr));
        y = zeros(length(sarr));
        for i = 1:length(sarr)
            x(i) = 2 * sqrt(b * log(lambdaMax / sarr(i)));
            y(i) = length(sarr) - i;
        end
        legendInfo{l} = ['L = ' num2str(L(l))];
        hold on
        scatter(x(:), y(:));
    end
    xlabel('2$\sqrt{bln(\frac{\lambda_{max}}{\lambda})}$', 'Interpreter', 'latex');
    ylabel('n($\lambda$)', 'Interpreter', 'latex');
    clx = (0:0.1:6);
    cly = besseli(0, (clx/2).^2);
    legendInfo{length(filenames) + 1} = ['CL'];
    hold on
    plot(clx, cly);    
    set(gca, 'YScale', 'log');
    xlabel('2$\sqrt{bln(\frac{\lambda_{max}}{\lambda})}$', 'Interpreter', 'latex');
    ylabel('n($\lambda$)', 'Interpreter', 'latex');
    legend(legendInfo);
    savefig('Fig3');
    hold off
end
