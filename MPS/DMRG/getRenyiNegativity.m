function rn = getRenyiNegativity(rhoT2, n, q)
    ind = find(rhoT2.Q{1} == q);
    if isempty(ind)
        rn = 0;
        return;
    end
    [~, v] = eig(rhoT2.data{ind});
    rn = sum(sum(v.^n));
end