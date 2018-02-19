function res = graph5()
    %
    Ls =     [100, 200, 400, 600, 800, 1000, 1500]; %, 2000];
    Delta = [0, -0.5, -0.9, 1];
    colors = ['r', 'b', 'g', 'k', 'y', 'c',  'm',  'r'];
    for l = 1 : length(Ls)
        for d = 1: length(Delta)
            s = load(strcat('spectrumN', int2str(Ls(l)), 'D', num2str(abs(Delta(d))), '_1024.mat'));
            marr = s.spectrum.keys();
            sz = zeros(1, length(marr));
            pm = zeros(1, length(marr));
            avgM = 0;
            avgM2 = 0;
            for m = 1:length(marr)
                sz(m) = str2num(marr{m}) / 2;
                vals = s.spectrum(marr{m});
                % keep only eigenvalues > 1e-9
                for i = 1:length(vals)
                    if vals(i) < 1e-9
                        vals = vals(1 : i-1);
                        break;
                    end
                    pm(m) = pm(m) + vals(i);
                end
                avgM = avgM + sz(m) * pm(m);
                avgM2 = avgM2 + sz(m)^2 * pm(m);
            end
            [sz, indices] = sort(sz);
            pm = pm(indices);
            f = fittype('b*exp(-(x^2)/(2*c))', 'independent', 'x', 'dependent', 'y');
    %         A = 0.01;
            B = pm(1);
            C = 2 * avgM2;
            fg = fit(sz.', pm.', f, 'StartPoint', [B C]);
            res.(strcat('L', int2str(Ls(l)), 'D', int2str(d), 'sz')) = sz;
            res.(strcat('L', int2str(Ls(l)), 'D', int2str(d), 'pm')) = pm;
            res.(strcat('L', int2str(Ls(l)), 'D', int2str(d), 'fg')) = fg;
            res.(strcat('L', int2str(Ls(l)), 'D', int2str(d), 'c2')) = C;
        end
    end
    % Plot fig 5
    for l = 1:length(Ls)
       h = plot(res.(strcat('L', int2str(Ls(l)), 'D3fg')), ...
                res.(strcat('L', int2str(Ls(l)), 'D3sz')), ...
                res.(strcat('L', int2str(Ls(l)), 'D3pm')));
       set([h(1) h(2)],'color',colors(l))
       legendVals(l) = h(2);
       legendInfo{l} = ['L = ' int2str(Ls(l))];
       hold on
    end
    set(gca, 'XScale', 'linear');
    set(gca, 'YScale', 'log');
    legend(legendVals, legendInfo);
    xlabel('m', 'Interpreter', 'latex');
    ylabel('$p_m(L/2)$', 'Interpreter', 'latex'); 
    savefig('Fig5');
    hold off;
    clearvars legendInfo legendVals;
    % Plot fig 5 inset
    x = Ls / 2;
    c2 = zeros(1, length(x));
    sigma = zeros(1, length(x));
    for l = 1 : length(Ls)
        c2(l) = res.(strcat('L', int2str(Ls(l)), 'D3c2'));
        sigma(l) = res.(strcat('L', int2str(Ls(l)), 'D3fg')).c;
    end
    scatter(x(:), sigma(:));
    hold on;
    scatter(x(:), c2(:));
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'linear');
    legend({'$\sigma^2$', 'C2'}, 'Interpreter','latex');
    xlabel('L/2', 'Interpreter', 'latex');
    savefig('Fig5Inset');
    hold off;
    % Plot fig 6(c)
    for d = 1:4
       h = plot(res.(strcat('L1500D', num2str(d), 'fg')), ...
                res.(strcat('L1500D', num2str(d), 'sz')), ...
                res.(strcat('L1500D', num2str(d), 'pm')));
       set([h(1) h(2)],'color',colors(d))
       legendVals(d) = h(2);
       legendInfo{d} = ['L = ' int2str(Ls(d))];
       hold on;
    end
    set(gca, 'XScale', 'linear');
    set(gca, 'YScale', 'log');
    legend(legendVals, legendInfo);
    xlabel('m', 'Interpreter', 'latex');
    ylabel('$p_m(L/2)$', 'Interpreter', 'latex'); 
    savefig('Fig5');
    hold off;
    clearvars legendInfo legendVals;
    % plot fig 6(a)
    for d = 1:4
        x = Ls / 2;
        c2 = zeros(1, length(x));
        sigma = zeros(1, length(x));
        for l = 1 : length(Ls)
            c2(l) = res.(strcat('L', int2str(Ls(l)), 'D', num2str(d), 'c2'));
            sigma(l) = res.(strcat('L', int2str(Ls(l)), 'D', num2str(d), 'fg')).c;
        end
        scatter(x(:), sigma(:), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors(d));
        legendInfo{d*2 - 1} = ['$\sigma^2$, $\Delta$ = ' num2str(Delta(d))];
        hold on;
        scatter(x(:), c2(:), 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors(d));
        legendInfo{d*2} = ['$C_2$, $\Delta$ = ' num2str(Delta(d))];
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'linear');
    legend(legendInfo, 'Interpreter','latex', 'Location','northwest');
    xlabel('L/2', 'Interpreter', 'latex');
    ylabel('$C_2$ and $\sigma^2$', 'Interpreter', 'latex');
    savefig('Fig6a');
    hold off;
end