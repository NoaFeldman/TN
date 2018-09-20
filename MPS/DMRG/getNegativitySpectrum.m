function spectrum = getNegativitySpectrum(E)
    % perform patial transpose and diagonalize:
    %      __
    %  a--|  |--c              __
    %     |  |       >>> a,d--|__|--c, b
    %  b--|__|--d
    E = contract(E, '23', getIdentity(E, 2, E, 3), '12');
    E = contract(E, '12', getIdentity(E, 1, E, 2), '12*');
    spectrum = [];
    for i = 1:length(E.data)
        [~, V] = eig(E.data{i});
        for i = 1 : length(V)
            v(i) = V(i, i);
        end
        spectrum = [spectrum v];
    end
end