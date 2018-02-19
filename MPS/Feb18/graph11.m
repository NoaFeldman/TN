function graph11()
    Ls =      [100, 200, 400, 600, 800, 1000, 1500, 2000];
    colors =  ['y', 'm', 'c', 'b', 'g', 'r',  'k',  'c'];
    symbols = ['>', 'v', '<', '^', 'd', 's',  'o',  '*'];
    Delta = [0, -0.5, -0.9, 1];
    for d = 1 : length(Delta)
        S1s = zeros(4, length(Ls));
        for l = 1 : length(Ls)
            if (d == length(Delta) & l == length(Ls))
                break;
            end
            s = load(strcat('spectrumN', int2str(Ls(l)), 'D', num2str(abs(Delta(d))), '_1024.mat'));
            marr = s.spectrum.keys();
            s1 = zeros(1, length(marr));
            for m = 1:length(marr)
                vals = s.spectrum(marr{m});
                % keep only eigenvalues > 1e-9
                for i = 1:length(vals)
                    if vals(i) < 1e-9
                        vals = vals(1 : i-1);
                        break;
                    end
                    s1(m) = s1(m) - vals(i) * log(vals(i));
                end
            end
        S1s(1:4, l) = s1(1:4);    
        end
        figure(1);
        scatter(Ls(:), S1s(1, :), 'markerEdgeColor', colors(d));
        hold on;
        legendInfo1{d} = ['$\Delta = $', num2str(Delta(d))];
        figure(2);
        B = 2;
        ceff = zeros(4, length(Ls));
        for i = 1:4
            for j = 1: length(Ls)
                ceff(i, j) = S1s(i, j) * 3 * B / log(Ls(j));
            end
        end
        scatter(Ls(:), ceff(1, :), 'markerEdgeColor', colors(d));
        hold on;
        legendInfo1{d} = ['$\Delta = $', num2str(Delta(d))];
    end
    figure(1);
    set(gca, 'XScale', 'log');
    legend(legendInfo1, 'Interpreter', 'latex');
    xlabel('L', 'Interpreter', 'latex');
    ylabel('$S^{(m=0)}_1$', 'Interpreter', 'latex'); 
    savefig('fig11aLeft');
    hold off;    
    figure(2);
    set(gca, 'XScale', 'log');
    legend(legendInfo1, 'Interpreter', 'latex');
    xlabel('L', 'Interpreter', 'latex');
    ylabel('$c_{eff}(m=0)$', 'Interpreter', 'latex'); 
    savefig('fig11aRight');
    hold off;
end
