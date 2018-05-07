function overlap = getOverlap(psia, psib, k)
    % Calculates < psia | psib >
    % k is the index of the site in mixed canonical form.
    if (~isempty(psib(1)))
        overlapLeft = getIdentity(psib(1), 1, [2 1]);
    else
        overlap = 0;
        return;
    end
    for i = 1:k
        overlapLeft = contract(contract(overlapLeft, 1, psib(i), 1), '12', psia(i), '12*');
    end
    if (~isempty(psib(length(psib))))
        overlapRight = getIdentity(psib(length(psib)), 3, [2 1]);
    else
        overlap = 0;
        return;
    end
    for i = length(psib):-1:k+1
        overlapRight =  contract(contract(overlapRight, 1, psib(i), 3), '13', psia(i), '32*');
    end
    overlap = getscalar(contract(overlapLeft, '12', overlapRight, '12'));
end