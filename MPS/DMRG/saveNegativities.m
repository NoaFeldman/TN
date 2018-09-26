function saveNegativities(psi, ls, ns, filename)
    L = length(psi);
    for i = 1:length(ns)
        n = ns(i);
        clear res;
        for j = 1:length(ls)
            l = ls(j);
            spectrum = containers.Map();
            [N1, N2] = getNegativityNs(psi, L/2 - l + 1, L/2, L/2 + l);
            rhoT2 = partiallyTransposedRDM(N1, N2);
            for q = 0:2:8
                spectrum(num2str(q)) = getRenyiNegativity(rhoT2, n, q);
            end
            res.(num2str(l)) = spectrum;
        end
        save(strcat(filename, '_', int2str(n)), res);
    end

end