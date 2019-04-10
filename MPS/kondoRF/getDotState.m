function res = getDotState(NRG, AV, dotVOp, dotCOp)
    chainStart = contract(contract(AV, 3, dotVOp, 2), '12', AV, '12*');
    chainStart.info.itags = {'*', ''};
    cLink = contract(contract(NRG(1).AK, 3, dotCOp, 2), '3', NRG(1).AK, '3*');
    currChain = contract(chainStart, '12', cLink, '13');
    for i = 2:length(NRG)-1
        currChain = contract(contract(currChain, 1, NRG(i).AK, 1), '13', NRG(i).AK, '13*');
    end
    currChain = contract(contract(currChain, 1, NRG(end).AT, 1), '13', NRG(end).AT, '13*');
    rho0 = getThermalState(NRG(end).HT, 0);             
    res = getscalar(contract(currChain, '12', rho0, '12'));  
end