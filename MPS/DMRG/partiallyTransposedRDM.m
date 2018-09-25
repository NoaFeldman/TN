function rhoT2 = partiallyTransposedRDM(N1, N2)
    % Fig 8 in https://arxiv.org/pdf/1605.00674.pdf
    % B1 (represented by E12) is empty.
    M = contract(N1, 2, N2, 1);
    % This action gets me out of memory - divide to sectors
    rho = contract(M, '13', M, '13*');
    rho2 = contract(contract(rho, 3, getIdentity(rho, 3, '-0', 'sA1'), '1', [1 2 4 3]), ...
        4, getIdentity(rho, 4, '-0', 'sA2'), 1);
    rhoT2 = contract(contract(rho, '14', getIdentity(rho, 1, rho, 4), '12*'), '12', ...
        getIdentity(rho, 2, rho, 3), '12*');    
%     rho = contract(contract(rho, 2, getIdentity(rho, 2, '-0', 'sA2'), '1*', [1 4 2 3]), ...
%         4, getIdentity(rho, 4, '-0', 'sA2'), 1);
%     rhoT2 = contract(contract(rho, '14', getIdentity(rho, 1, rho, 4), '12*'), '12', ...
%         getIdentity(rho, 2, rho, 3), '12*');
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