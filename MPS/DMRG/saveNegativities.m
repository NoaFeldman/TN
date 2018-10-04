function saveNegativities(psi, ls, ns, filename)
    maxNumCompThreads(4);
    L = length(psi);
    for j = 1:length(ls)
        l = ls(j);
        spectrum = containers.Map();
        [N1, N2] = getNegativityNs(psi, L/2 - l + 1, L/2, L/2 + l); % 5*L/8); %
        rhoT2 = partiallyTransposedRDM(N1, N2);
        for q = 0:2:12
            spectrum(num2str(q)) = getRenyiNegativity(rhoT2, ns, q);
        end
        res.(strcat('l',num2str(l))) = spectrum;
    end
    save(filename, 'res', 'ns', 'ls');
end
