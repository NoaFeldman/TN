function negs = arrangeNegativityDMRG(filenames, ls)
    qs = -6:6;
    rn3 = zeros(length(ls), length(qs));
    for i = 1:length(ls)
        l = ls(i);
        specMap = getSpecMap(filenames, l);
        for j = 1:length(qs)
            q = qs(j);
            spectrum = specMap(int2str(2*abs(q)));
            rn3(i, j) = sum(spectrum.^3);
        end
        normalized = rn3(i, :) ./ sum(rn3(i, :));
        f = fittype('1/sqrt(2*pi*c)*exp(-x^2/(2*c))', 'independent', 'x', 'dependent', 'y');
        [fg, gof] = fit(qs.', normalized.', f, 'StartPoint', [0.1]);
        plot(fg, qs.', normalized); pause(0.1);
        sig2(i) = fg.c;
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
    
