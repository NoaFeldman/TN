function arranged = arrangeNegativityDMRGRenyi(map, n)
    keys = map.keys;
    vals = map.values;
    for i = 1 : length(keys)
        qs(i) = str2num(keys{i})/2;
        rnq(i) = sum(vals{i}.^n);
    end
    [qs, sortInds] = sort(qs);
    rnq = rnq(sortInds);
    alphas = -3.14:0.01:3.14;
    rnalpha = zeros(1, length(alphas));
    for i = 1:length(qs)
        rnalpha = rnalpha + exp(1i * qs(i) .* alphas).* rnq(i);
        if qs(i) ~= 0
            rnalpha = rnalpha + exp(1i * (-qs(i)) .* alphas).* rnq(i);
        end
    end
    arranged.qs = qs;
    arranged.rnq = rnq;
    arranged.alphas = alphas;
    arranged.rnalpha = rnalpha;
end
