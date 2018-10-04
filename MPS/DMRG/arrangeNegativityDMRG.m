function negs = arrangeNegativityDMRG(filename)
    q = -6:6;
    orig = load(filename);
    f = fittype('b*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
    for i = 1:length(orig.ns)
        n = orig.ns(i);
        sigma = zeros(1, length(orig.ls));
        rFull = zeros(1, length(orig.ls));
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
            [fg, gof] = fit(q.', y.', f, 'StartPoint', [rn 1]);
            plot(fg, q, y); pause(0.1);
            sigma(j) = fg.c;
            rFull(j) = fg.b * sqrt(2 * pi * fg.c); 
        end
    end
end