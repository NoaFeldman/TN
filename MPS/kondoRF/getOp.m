function op = getOp(NRG, k, op0)
    kp = 1;
    op = op0;
    while kp <= k
        op = getOpStep(NRG, kp, op);
        kp = kp + 1;
    end
end


function op = getOpStep(NRG, k, opKMin1)
    if k == 1
        op = contract(NRG(k).AK, '13*', contract(NRG(k).AK, 3, opKMin1, 2), '13'); 
    else
        op = contract(NRG(k).AK, '13*', contract(opKMin1, 2, NRG(k).AK, 1), '13');
    end
end