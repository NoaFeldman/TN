function rn = getRenyiNegativity(rhoT2, ns, q)
    ind = find(rhoT2.Q{1} == q);
    if isempty(ind)
        rn = 0;
        return;
    end
    [~, v] = eig(rhoT2.data{ind});
    for i = 1:length(ns)
        n = ns(i);
        rn(i) = sum(sum(v.^n));
    end
end