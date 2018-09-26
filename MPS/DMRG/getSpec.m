function spec = getSpec(m)
    % Get spectrum of a 2-dim QSpace
    spec = [];
    for i = 1:length(m.data)
        [~, v] = eig(m.data{i});
        for j = 1:length(v)
            spec = [spec v(j, j)];
        end
    end
end