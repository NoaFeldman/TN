function res = findSEps(t, s, epsP)
    NA = -2:0;
    epsS = 1:1000;
    chi = zeros(1, length(eps));
    for na = 1:length(NA)
        sTheo = stheo(t.', NA(na), epsP, reshape(epsS,1,1,[]));
        sMat = repmat(s(na, :).', 1, length(epsP), length(epsS));
        chi = chi + (sum((sMat - sTheo).^2, 1));
        res.(strcat('stheo', int2str(na))) = sTheo;
        res.(strcat('chi', int2str(na))) = chi;
    end
end