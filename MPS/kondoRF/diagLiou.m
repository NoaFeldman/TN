function [rAlphas, lAlphas, lambdaAlphas] = diagLiou(Liou, rhoT)
    calcOverlap = @(rho1, rho2) trace(rho1'*rho2);
    applyLiou = @(rho, Liou) applyLiouvillian(rho, Liou);
    newStateArr = @(n) QSpace(n);
       
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
    