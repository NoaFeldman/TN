function graph9()
    Ls =      [100, 200, 400, 600, 800, 1000, 1500, 2000];
    colors =  ['y', 'm', 'c', 'b', 'g', 'r',  'k',  'c'];
    symbols = ['>', 'v', '<', '^', 'd', 's',  'o',  '*'];
    Delta = [0, -0.5, -0.9, 1];
    for d = 1 : length(Delta)
        for l = 1 : length(Ls)
            if (d == length(Delta) & l == length(Ls))
                break;
            end
            s = load(strcat('spectrumN', int2str(Ls(l)), 'D', num2str(abs(Delta(d))), '_1024.mat'));
            marr = s.spectrum.keys();
            ms = 0;
            wrapper = 0;
            for m = 1:length(marr)
                vals = s.spectrum(marr{m});
                if (marr{m} == '0')
                    lambdaMax = vals(1);
                end
                % keep only eigenvalues > 1e-9
                for i = 1:length(vals)
                    if vals(i) < 1e-9
                        vals = vals(1 : i-1);
                        break;
                    end
                    if (i == 1)
                        index = 1;
                        if (m > 1)
                            index = length(ms) + 1;
                        end
                        ms(index) = str2num(marr{m}) / 2;
                        wrapper(index) = -log(vals(i) / lambdaMax);
                    end
%                     if (d == 2)
%                         x(i) = str2num(marr{m}) / 2;
%                         y(i) = -log(vals(i) / lambdaMax);
%                      figure(1)
%                      scatter(x, y, symbols(l), 'MarkerEdgeColor', colors(l));
%                      hold on;
%                     end
                end
            end
            f = fittype('b*x^2', 'independent', 'x', 'dependent', 'y');
            fp = fit(ms.', wrapper.', f, 'StartPoint', [pi^2]);
            b(l) = fp.b;
%             if (d == 2)
%                 figure(1)
%                 p = plot(fp, ms, wrapper);
%                 hold on;
%                 legendVals1(l) = [p(2)];
%                 legendInfo1{l} = ['L = ' int2str(Ls(l))];
%                 set([p(1) p(2)],'color',colors(l));
%             end
            clearvars ms wrapper
        end
        figure(2);
        scatter(Ls(:), b(:), 'markerEdgeColor', colors(d));
        legendInfo2{d} = ['\Delta = ' num2str(Delta(d))];
        hold on;
    end
%     figure(1);
%     legend(legendVals1, legendInfo1);
%     xlabel('m', 'Interpreter', 'latex');
%     ylabel('$-ln(\lambda_i/\lambda_max)$', 'Interpreter', 'latex'); 
%     savefig('figs/fig9');
%     hold off;
    figure(2);
    set(gca, 'XScale', 'log');
    legend(legendInfo2);
    xlabel('L', 'Interpreter', 'latex');
    ylabel('$\frac{m^2}{S^{(m)}_\infty - S_\infty}$', 'Interpreter', 'latex'); 
    savefig('figs/fig10');
    hold off;
end
            
            