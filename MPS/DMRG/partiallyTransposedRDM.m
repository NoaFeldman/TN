function rhoT2 = partiallyTransposedRDM(N1, N2)
    % Fig 8 in https://arxiv.org/pdf/1605.00674.pdf
    % B1 (represented by E12) is empty.
    M = contract(N1, 2, N2, 1);
    % This action gets me out of memory - divide to sectors
    rho = contract(M, '13', M, '13*');
    clear M;
    % For symmetric cases (h = 0), don't save matrices corresponding to negative DeltaN.
    inds = find(rho.Q{1} - rho.Q{4} >= 0);
    for i = 1:4
        rho.Q{i} = rho.Q{i}(inds);
    end
    rho.data = rho.data(inds);
    
    rhoT2 = rho;
    rhoT2.Q{2} = -rho.Q{4};
    rhoT2.Q{4} = -rho.Q{2};
    for i = 1:length(rho.data)
        rhoT2.data{i} = permute(rho.data{i}, [1 4 3 2]);
    end
    clear rho;
    % getIdentity doesn't work on complex QSpaces for some reason, so I use
    % temp for creating the index unifying QSpaces.
    temp = rhoT2;
    for i = 1:length(temp.data)
        temp.data{i} = real(temp.data{i});
    end
    rhoT2 = contract(contract(rhoT2, '12', getIdentity(temp, 1, temp, 2), '12*'), '12',...
        getIdentity(temp, 3, temp, 4), '12');
end