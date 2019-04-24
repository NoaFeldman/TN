function [LiouMatD, basisD, rhoVecD] = getDenseBasis(rhoVec, LiouMat)
    blocks = findMatrixBlocks(LiouMat);
    rhoInds = find(rhoVec ~= 0);
    bInds = [];
    for b = 1:length(blocks)
        if length(union(blocks{b}, find(rhoVec ~= 0))) < length(blocks{b}) + length(find(rhoVec ~= 0))
            bInds(end+1) = b;
        end
    end
    basisD = [];
    for i = 1:length(bInds)
        basisD = [basisD; blocks{bInds(i)}];
    end
    LiouMatD = LiouMat(basisD, basisD);
    rhoVecD = rhoVec(basisD);
end