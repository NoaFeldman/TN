function [rAlphas, lAlphas, lambdaAlphas, rhoSS, Liou, rhoT] = diagLiou(H, Sminus, SplusSminus, EScale, gamma, T)
    rhoT = getThermalState(H, 1);
    RI.H = H;
    RI.Sminus = Sminus;
    RI.SplusSminus = SplusSminus;
    rhoSS = RabiKondo.getSteadyState(RI, gamma, EScale);
    calcOverlap = @(rho1, rho2) trace(rho1'*rho2);
    applyLiou = @(rho, Liou) applyLiouvillian(rho, Liou);
    newStateArr = @(n) QSpace(n);
    rhoT = rhoT / sqrt(calcOverlap(rhoT, rhoT));
    
    Liou = getLiouvillian(H, gamma/EScale, Sminus, SplusSminus);
    
    if 1
        NKeep = 1024;

        [~, ~, lambdaAlphas, rAlphas, lAlphas] = ...
            arnoldi(rhoT, Liou, NKeep, applyLiou, calcOverlap, newStateArr);        
    else
        NKeep = 116;

        [~, ~, lambdaAlphas, rAlphas, lAlphas] = ...
            arnoldi(rhoT, Liou, NKeep, applyLiou, calcOverlap, newStateArr, 'brut');
    end
end
    