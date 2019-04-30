function [rAlphas, lDagAlphas, lambdaAlphas, rhoSS] = diagVecedLiou(vecedLiou, id, idIn)
    % Add for each block in vecedLiou
    currLiou = contract(vecedLiou, '12', id, '12*');
    [evecsR, evals, evecsL] = eig(currLiou.data{1});
    lambdaAlphas = diag(evals);
    
    rAlphas = currLiou;
    rAlphas.data = {evecsR};
    [~, ind] = min(lambdaAlphas);
    rhoSS = rAlphas;
    rhoSS.data = {evecsR(:, ind)'};
    rAlphas.info.itags{1} = 'alpha';
    rAlphas = contract(rAlphas, 2, id, 3);
    rAlphas = contract(rAlphas, 3, idIn', '1');
    rhoSS.info.itags{1} = 'alpha';
    rhoSS = contract(rhoSS, 2, id, 3);
    rhoSS.Q = rhoSS.Q(2:3);
    rhoSS.info.itags = rhoSS.info.itags(2:3);
    rhoSS.data{1} = reshape(rhoSS.data{1}, [size(rhoSS.data{1}, 2), size(rhoSS.data{1}, 3)]);
    rhoSS = contract(rhoSS, 2, idIn', '1');
    lDagAlphas = currLiou;
    lDagAlphas.data = {conj(evecsL)};
    lDagAlphas.info.itags{1} = 'alpha';
    lDagAlphas = contract(lDagAlphas, 2, id, 3, [1 3 2]);
    lDagAlphas = contract(lDagAlphas, 3, idIn', '1');    
end