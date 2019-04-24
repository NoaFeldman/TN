function op = getOp(NRG, op0)
    k = 1;
    op = op0;
    while k <= length(NRG) - 1
        if k == 1
            op(k) = contract(NRG(k).AK, '13*', contract(NRG(k).AK, 3, op0, 2), '13');
        else
            op(k) = contract(NRG(k).AK, '13*', contract(op(k-1), 2, NRG(k).AK, 1), '13');
        end
        k = k + 1;
    end
end