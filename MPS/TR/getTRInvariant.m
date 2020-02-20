function z2 = getTRInvariant(psi, A1Start, A1End, A2End)
    k = length(psi);
    while k > A1End
        [psi, k, truncErr] = shiftWorkingSite(psi, k, '<<', {});
    end
    [UT, utMatrix] = getTimeReversealUnitary();
    
    ut1 = UT;
    ut1.info.itags = {strcat(int2str(A1Start), 's'), strcat(int2str(A1Start), 's*')};
    A1Curr = contract(contract(psi(A1Start), 2, ut1, 2), 3, psi(A1Start), '2*');
    for k = A1Start + 1 : A1End
        utk = UT;
        utk.info.itags = {strcat(int2str(k), 's'), strcat(int2str(k), 's*')};
        temp = contract(psi(k), 2, utk, 2);
        A1Curr = contract(contract(A1Curr, 2, temp, 1), '35', psi(k), '12*', [1, 3, 2, 4]);
    end
    A1Full = contract(A1Curr, '13', A1Curr, '31');
    
    A2Curr = contract(psi(A1End + 1), '2', psi(A1End + 1), '2*');
    for k = A1End + 2 : A2End
        A2Curr = contract(contract(A2Curr, 2, psi(k), 1), '34', psi(k), '12*', [1, 3, 2, 4]);
    end
    A2Full = contract(A2Curr, '24', A2Curr, '42');
%     A2Full.Q{2} = -A2Full.Q{2};
%     A2Full.Q{3} = -A2Full.Q{3};
    A2Full.info.itags{3} = A2Full.info.itags{2};
    A2Full.info.itags{2} = strcat(A2Full.info.itags{2}, '*');
    z2 = getscalar(contract(A1Full, '1234', A2Full, '1324'));
end
            