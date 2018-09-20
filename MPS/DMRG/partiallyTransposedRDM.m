function rhoT2 = partiallyTransposedRDM(N1, N2)
    % Fig 8 in https://arxiv.org/pdf/1605.00674.pdf
    % B1 (represented by E12) is empty.
    M = contract(N1, 2, N2, 1);
    % This action gets me out of memory - divide to sectors
    rho = contract(M, '13', M, '13*');
    rho.info.itags(3:4) = {'sA1Prime', 'sA2Prime'};
    rho.Q{3} = -rho.Q{3};
    rho.Q{4} = -rho.Q{4};
    rhoT2 = contract(contract(rho, '14', getIdentity(rho, 1, rho, 4), '12*'), '12', ...
        getIdentity(rho, 2, rho, 3), '12*');
end