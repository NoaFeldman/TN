function E = getNegativityE(psi, startIndex, endIndex)
    % compute E matrices from figure 6 at https://arxiv.org/pdf/1605.00674.pdf
    % Working site is expected to be psi(endIndex).
    % startIndex and endIndex start from the right.
    ERight = contract(QSpace(psi(endIndex)), '2', QSpace(psi(endIndex)), '2*');
    for currIndex = endIndex - 1 : -1 : startIndex
        %  --O-- 1 --|--
        %    |       |  
        %    2       | 
        %    |       |
        %  --O-- 2 --|--
        temp = contract(QSpace(psi(currIndex)), 3, ERight, 1);
        ERight = contract(temp, '24', psi(currIndex), '23*', [1 2 4 3]);
    end
    E = ERight;
end