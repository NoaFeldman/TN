function res = getDotState(NRG, AV, dotVOp, dotCOp)
    chainStart = contract(contract(AV, 2, dotVOp, 2), 2, AV, '2*');
    chainStart.info.itags = {'*', ''};
    cLink = contract(contract(NRG(1).AK, 3, dotCOp, 2), '3', NRG(1).AK, '3*');
    currChain = contract(chainStart, '12', cLink, '13');
    for i = 2:length(NRG)-1
        currChain = contract(contract(currChain, 1, NRG(i).AK, 1), '13', NRG(i).AK, '13*');
    end
    if (isempty(currChain))
        res = 0;
    else
        currChain = contract(currChain, '12', getIdentity(currChain, 1), '12');
        res = getscalar(currChain);
    end   
end