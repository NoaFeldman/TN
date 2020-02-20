function saveRQ34(inDirName, outDirName, ti, tf, l)
    for t = ti:2:tf
        f = load(strcat(inDirName, '/psiAtStep', int2str(t)));
        opts = {'Nkeep', 512};
        L = length(f.psi);
        [N1, N2, truncErrs] = getNegativityNs(f.psi, L/2 - l + 1, L/2, L/2 + l, opts);
        
        spec = containers.Map();
        rhoT2Q2 = partiallyTransposedRDM(N1, N2, 2);
        [~, V] = eig(rhoT2Q2.data{1});
        spec('2') = diag(V);
        rhoT2Q3 = partiallyTransposedRDM(N1, N2, 3);
        [~, V] = eig(rhoT2Q3.data{1});
        spec('3') = diag(V);
        rhoT2Q4 = partiallyTransposedRDM(N1, N2, 4);
        [~, V] = eig(rhoT2Q4.data{1});
        spec('4') = diag(V);
        save(strcat(outDirName, '/step', int2str(t)), 'spec');
    end
end
        
