function saveNegativitySpectrum(psi, ls, mode, filename)
    L = length(psi);
    for j = 1:length(ls)
        l = ls(j);
        renyis = containers.Map();
        if strcmp(mode, 'symm')
            v2 = L/2 + l;
        elseif strcmp(mode, 'asymm')
            v2 =  5*L/8;
        end
        opts = {'Nkeep', 128};
        [N1, N2, truncErrs] = getNegativityNs(psi, L/2 - l + 1, L/2, v2, opts);
        rhoT2 = partiallyTransposedRDM(N1, N2);
        for q = 0:2:18
            renyis(num2str(q)) = getSpectrum(rhoT2, q);
        end
        res.(strcat('l',num2str(l))) = renyis;
        res.truncErrs = truncErrs;
    end
    save(filename, 'res', 'ls');
end

function spec = getSpectrum(rhoT2, q)
    ind = find(rhoT2.Q{1} == q);
    if isempty(ind)
        spec = [];
    else
        [~, V] = eig(rhoT2.data{ind});
        spec = diag(V);
    end
end
