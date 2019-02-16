function rhoT2 = partiallyTransposedRDM(N1, N2, q)
    % Fig 8 in https://arxiv.org/pdf/1605.00674.pdf
    % B1 (represented by E12) is empty.
    M = contract(N1, 2, N2, 1);
    % This action gets me out of memory - divide to sectors
    chargeA1 = getIdentity(M, 2);
    chargeA2 = getIdentity(M, 4);
    rhoQ = QSpace();
    for i = 1:length(chargeA1.Q{1})
        j = find(-chargeA2.Q{2} == q*2 - chargeA1.Q{1}(i));
        proj1 = getProjector(chargeA1, i);
        proj2 = getProjector(chargeA2, j);
        MTemp = contract(contract(proj1, 1, M, 2, [2 1 3 4]), 4, ...
                         proj2, 2);
        rhoQ = rhoQ + contract(MTemp, '13', MTemp, '13*');
    end
    
%     rho = contract(M, '13', M, '13*');
    clear M;
    % For symmetric cases (h = 0), don't save matrices corresponding to negative DeltaN.
%     inds = find(rho.Q{1} - rho.Q{4} == 4);
%     inds = [inds; find(rho.Q{1} - rho.Q{4} == 6)];
%     for i = 1:4
%         rho.Q{i} = rho.Q{i}(inds);
%     end
%     rho.data = rho.data(inds);
     
    rhoT2 = rhoQ;
    rhoT2.Q{2} = -rhoQ.Q{4};
    rhoT2.Q{4} = -rhoQ.Q{2};
    for i = 1:length(rhoQ.data)
        rhoT2.data{i} = permute(rhoQ.data{i}, [1 4 3 2]);
    end
    clear rho;
    rhoT2 = contract(contract(rhoT2, '12', getIdentity(real(rhoT2), 1, real(rhoT2), 2), '12*'), '12',...
        getIdentity(real(rhoT2), 3, real(rhoT2), 4), '12');
end

function proj = getProjector(chargeA, i)
    proj = chargeA;
    proj.data = proj.data(i);
    proj.Q{1} = proj.Q{1}(i);
    proj.Q{2} = proj.Q{2}(i);
end