function [rAlphas, lAlphas, lambdaAplhas] = diagLiou(H, Sminus, SplusSminus, EScale, gamma, T)
    rhoT = getThermalState(H, T);
    RI.H = H;
    RI.Sminus = Sminus;
    RI.SplusSminus = SplusSminus;
    rhoSS = RabiKondo.getSteadyState(RI, gamma, EScale);
    calcOverlap = @(rho1, rho2) trace(rho1'*rho2);
    applyLiou = @(rho, Liou) applyLiouvillian(rho, Liou);
    newStateArr = @(n) QSpace(n);
    
    Liou = getLiouvillian(H, gamma, Sminus, SplusSminus);
    
    NKeep = 1024;
    
    [~, ~, lambdaAplhas, rAlphas, lAlphas] = ...
        arnoldi(rhoT, Liou, NKeep, applyLiou, calcOverlap, newStateArr);
end
    