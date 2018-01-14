function HPairs = getHPairs(H, dtReal, dtIm)
    % TODO
    N = length(H.single);
    HPairs = QSpace(N - 1);
    for k = 1 : N - 1
        HPairs(k) = contract(H.l2r(k), 3, H.r2l(k+1), 3);
        pairID = getIdentity(H.single(k), 1, H.identity(k+1), 1);
        pairID = contract(pairID, 3, pairID, '3*', [3 1 4 2]);
        HPairs(k) = HPairs(k) + contract(pairID, 2, H.single(k), 1, [1 4 2 3]);
        if (k == N - 1)
            HPairs(k) = HPairs(k) + contract(pairID, 4, H.single(k+1), 1);
        end
        %exponentiate
        for l=1:length(HPairs(k).data)
            HPairs(k).data{l} = expm(-0.5 * complex(0, 1) * complex(dtReal, dtIm).* HPairs(k).data{l});
        end
    end
end
