function [rAlphas, lDagAlphas, lambdaAlphas, rhoSS] = diagVecedLiou(vecedLiou, id)
    % Add for each block in vecedLiou
    currLiou = contract(vecedLiou, '12', id, '12*')';
    [evecsR, evals, evecsL] = eig(currLiou.data{1});
    lambdaAlphas = diag(evals);
    
    rAlphas = currLiou;
    rAlphas.data = {evecsR};
    [~, ind] = min(lambdaAlphas);
    rhoSS = rAlphas;
    rhoSS.data = {evecsR(:, ind)};
    rAlphas.info.itags{2} = 'alpha';
    rAlphas = contract(rAlphas, '1*', id, '3');
    rhoSS.info.itags{2} = 'alpha';
    rhoSS = contract(rhoSS, '1*', id, '3');
    rhoSS.Q = rhoSS.Q(2:3);
    rhoSS.info.itags = rhoSS.info.itags(2:3);
    rhoSS.data{1} = reshape(rhoSS.data{1}, [size(rhoSS.data{1}, 2), size(rhoSS.data{1}, 3)]);
    
    n = length(evecsL);
    evecsLDag = zeros(n);
    for i = 1:n
        evecsLDag(:, i) = reshape(reshape(evecsL(:, i), [sqrt(n) sqrt(n)])', [n 1]);
    end
    lDagAlphas = currLiou;
    lDagAlphas.data = {evecsLDag};
    lDagAlphas.info.itags{2} = 'alpha';
    lDagAlphas = contract(lDagAlphas, '1*', id, '3');
end