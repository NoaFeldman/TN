function spec = processRDMSpectrum(filename)
% TODO unify with arrangeDMRGResults
    sz = -10:2:10;
    spec.s = zeros(length(sz), 1); 
    spec.p = zeros(length(sz), 1);
    spec.sFull = 0;
    file = load(filename);
    k = keys(file.spectrum);
    val = values(file.spectrum);
    for i = 1:length(val)
        currSz = str2num(k{i});
        ind = find(sz == currSz);
        pNA = sum(val{i});
        dn = 1e-2;
        sNA = -(sum(val{i}.^(1+dn)) - sum(val{i}.^(1-dn))) / (2 * dn);
        if (~isempty(ind))
            spec.s(ind) = sNA;
            spec.p(ind) = pNA;
        end
        spec.sFull = spec.sFull + sNA;
    end
    spec.alphas = -3.14:0.01:3.14;
    for i = 1:length(spec.alphas)
        spec.salpha(i) = sum(spec.p.*exp(1i .* ((sz./2).') .* spec.alphas(i)));
    end
end