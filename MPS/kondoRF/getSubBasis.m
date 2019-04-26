function [LiouMatSub, basisSub, evecsRSub, evecsLSub, evalsSub] = getSubBasis(rhoVec, LiouMat)
    blocks = findMatrixBlocks(LiouMat);
    bInds = [];
    for b = 1:length(blocks)
        if length(union(blocks{b}, find(rhoVec ~= 0))) < length(blocks{b}) + length(find(rhoVec ~= 0))
            bInds(end+1) = b;
        end
    end
    basisSub = [];
    for i = 1:length(bInds)
        basisSub = [basisSub; blocks{bInds(i)}];
    end
    LiouMatSub = LiouMat(basisSub, basisSub);
    [evecsRTmp, evalsTmp, evecsLTmp] = eig(LiouMatSub);
    evalsSub = diag(evalsTmp);

    evecsRSub = zeros(length(LiouMat), length(basisSub));
    evecsRSub(basisSub, :) = evecsRTmp;
    evecsLSub = zeros(length(LiouMat), length(basisSub));
    evecsLSub(basisSub, :) = evecsLTmp;
end