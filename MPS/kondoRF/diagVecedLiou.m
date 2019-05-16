function [rAlphas, lDagAlphas, lambdaAlphas, rhoSS] = diagVecedLiou(vecedLiou, id)
    % get a matrix:
    %        |
    %   ____id_____
    %  |           |
    %  |           | 
    %  |_vecedLiou_| 
    %        |     
    currLiou = contract(vecedLiou, '12', id, '12*');
    if length(currLiou.data) > 1
        disp('diagVecedLiou: vecedLiou is not a single block???');
    end
    [evecsR, evals, evecsL] = eig(currLiou.data{1});
    evecsL = normalizeEvecs(evecsR, evecsL);
    lambdaAlphas = diag(evals);
    
    % 1 rank tensor is not right for QSpace, so we crate a 2-rank tensor
    % with one index of the vector and one alpha index (alpha counts the
    % eigenvectors).
    rAlphas = newQSpace({[0 0], [0 0]}, {evecsR}, {'rho*', 'alpha*'});
    rAlphas = contract(id, '3*', rAlphas, 1);
    lDagAlphas = newQSpace({[0 0], [0 0]}, {evecsL}, {'rho*', 'alpha*'});
    lDagAlphas = contract(id, '3*', lDagAlphas, 1);
    [~, ind] = min(lambdaAlphas);
    rhoSS = newQSpace({[0 0], [0 0]}, {evecsR(:, ind)}, {'rho*', 'alpha*'});
    rhoSS = contract(id, '3*', rhoSS, 1);
    % Lose the alpha leg
    rhoSS.info.itags = rhoSS.info.itags(1:2);
    rhoSS.Q = rhoSS.Q(1:2);
end

function evecsL = normalizeEvecs(evecsR, evecsL)
    for i = 1:length(evecsR)
        rMat = reshape(evecsR(:, i), [sqrt(length(evecsR)) sqrt(length(evecsR))]);        
        lDagMat = reshape(evecsL(:, i), [sqrt(length(evecsR)) sqrt(length(evecsR))])';
        evecsL(:, i) = evecsL(:, i) / trace(rMat * lDagMat);
    end
end