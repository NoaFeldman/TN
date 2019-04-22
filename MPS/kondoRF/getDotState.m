function res = getDotState(NRGPath, AV, AC, dotVOp, dotCOp, T)
    [~, ~, ~, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    distr = getDistribution(NRGPath, AV, AC, dotVOp, dotCOp, IS.E, T, 'overlap');
    res = distr(end);
end