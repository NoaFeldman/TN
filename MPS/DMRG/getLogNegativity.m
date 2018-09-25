function epsilon = getLogNegativity(rhoT2)
    epsilon = log(sum(abs(getSpec(rhoT2))));
end

function spec = getSpec(m)
    spec = [];
    for i = 1:length(m.data)
        [~, v] = eig(m.data{i});
        for j = 1:length(v)
            spec = [spec v(j, j)];
        end
    end
end