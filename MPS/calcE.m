HPair = contract(S(2), 3, S(2)', 3)  + contract(S(3), 3, S(3)', 3);
SZ2 = S(1);
SZ2.info.itags{3} = '';
HPair = HPair + JZ * contract(S(1), 3, SZ2, 3); 
opts = {'NKeep', 1024};
E = 0;
for i = 1:length(psi)-1
    curr = psi;
    k = length(curr);
    while k > i+1
        [curr, k] = shiftWorkingSite(curr, k, '<<');
    end
    M = contract(curr(i), 3, curr(i+1), 1);
    M = contract(M, '23', HPair, '24', [1 3 4 2]);
    [r, l, I] = orthoQS(M, [1 2], '>>');
    curr(i) = QSpace(r);
    curr(i+1) = QSpace(l);
    while k < length(curr)
        [curr, k] = shiftWorkingSite(curr, k, '>>');
    end
    E = E + getOverlap(psi, curr, length(psi));
end