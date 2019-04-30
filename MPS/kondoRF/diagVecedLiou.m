function [rAlphas, lDagAlphas, lambdaAlphas] = diagVecedLiou(vecedLiou, id)
    [evecsR, evals, evecsL] = eig(vecedLiou.data{1});
    lambdaAlphas = diag(evals);
    
    rAlphas = vecedLiou;
    rAlphas.data = {evecsR};
    rAlphas.info.itags{1} = 'alpha';
    rAlphas = contract(rAlphas, 2, id, 3);
    lDagAlphas = vecedLiou;
    lDagAlphas.data = {conj(evecsL)};
    lDagAlphas.info.itags{1} = 'alpha';
    lDagAlphas = contract(lDagAlphas, 2, id, 3, [1 3 2]);
    
end