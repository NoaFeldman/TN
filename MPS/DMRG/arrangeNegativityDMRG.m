function negs = arrangeNegativityDMRG(filename)
    q = -6:6;
    orig = load(filename);
    for i = 1:length(orig.ns)
        n = orig.ns(i);
        sigma = zeros(1, length(orig.ls));
        for j = 1:length(orig.ls)
            l = orig.ls(j);
            y = zeros(1, length(q));
            map = orig.res.(strcat('l', int2str(l)));
            for k = 1:length(q)
                vals = map(int2str(2*(q(k))));
                if vals == 0
                    y(k) = 0;
                else
                    y(k) = vals(i);
                end
            end
            rn = sum(y);
            f = fittype(strcat(num2str(rn), '*1/sqrt(2*pi*c)*exp(-x^2/(2*c))'), 'independent', 'x', 'dependent', 'y');
            [fg, gof] = fit(q.', y.', f, 'StartPoint', [0.1]);
            plot(fg, q, y); pause(0.1);
            sigma(j) = fg.c;
        end
    end
end