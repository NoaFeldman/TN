function omegas = testingOmegas(LiouB, rhoTB, SminusB, SplusB)
    [LiouMat, rhoVec] = vecing(LiouB, rhoTB);
    [evecsR, evals, evecsL] = eig(LiouMat);
    evals = diag(evals);
    rAlphas = {};
    lDagAlphas = {};
    for i = 1:length(evals)
        rAlphas{i} = reshape(evecsR(:, i), [sqrt(length(evecsR(:, i))) sqrt(length(evecsR(:, i)))]);
        lDagAlphas{i} =  reshape(evecsL(:, i), [sqrt(length(evecsR(:, i))) sqrt(length(evecsR(:, i)))])';
        lDagAlphas{i} = lDagAlphas{i}/trace(rAlphas{i} * lDagAlphas{i});
    end
    [~, ind] = min(evals);
    rhoSS = reshape(evecsR(:, ind), [sqrt(length(evecsR(:, i))) sqrt(length(evecsR(:, i)))]);
    SminusMat = SminusB.data{1};
    SplusMat = SplusB.data{1};
    
    for alpha = 1:length(evals)
        omegas(alpha) = trace(SminusMat * rAlphas{alpha}) * ...
                        trace(lDagAlphas{alpha} * (SplusMat * rhoSS - rhoSS * SplusMat));
    end
end