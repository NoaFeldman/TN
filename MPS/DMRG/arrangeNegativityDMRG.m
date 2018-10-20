function negs = arrangeNegativityDMRG(filenames, ls)
    qs = -6:6;
    rn3 = zeros(length(ls), length(qs));
    ns = 1:6;
    lRight = 16;
    for i = 1:length(ns);
        n = ns(i);
        for j = 1:length(ls)
            l = ls(j);
            specMap = getSpecMap(filenames, l);
            for k = 1:length(qs)
                q = qs(k);
                spectrum = specMap(int2str(2*abs(q)));
                rn(j, k) = sum(spectrum.^n);
            end
            normalized = rn(j, :) ./ sum(rn(j, :));
            f = fittype('1/sqrt(2*pi*c)*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
            [fg, gof] = fit(qs.', normalized.', f, 'StartPoint', [0.1]);
            plot(fg, qs.', normalized); pause(0.5);
            sig2(j) = fg.c;
        end
        scatter(ls, sig2);
        xlabel('$l_1$', 'Interpreter', 'latex');
        ylabel('$Var_3(Q)$', 'Interpreter', 'latex');
        title(strcat('Variance of $R_n(Q)/R_n$, $n = ', int2str(n), '$, $l_2 = 16$'), 'Interpreter', 'latex');
        saveas(gcf, strcat('t0L128n',int2str(n), '.svg'));
        scatter(log(ls.^2./(ls + 16)), sig2);
        scatter(log(ls(4:12).^2./(ls(4:12) + 16)), sig2(4:12));
        hl = lsline;
        B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
        B * n * pi^2;
    end
end

function specMap = getSpecMap(filenames, l)
    for k = 1:length(filenames)
        load(filenames{k});
        if (isfield(res, strcat('l', int2str(l))))
            specMap = res.(strcat('l', int2str(l)));
        end
    end
%     specMap = containers.map();
end
    
