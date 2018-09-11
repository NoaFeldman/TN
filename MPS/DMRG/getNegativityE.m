function E = getNegativityE(psi, startIndex, endIndex)
    % compute E matrices from figure 6 at https://arxiv.org/pdf/1605.00674.pdf
    % psi is expected to be left canonical, B is the middle of the chain.
    % startIndex and endIndex start from the right.
    if (endIndex == length(psi))
        ERight = contract(psi(endIndex), '23', psi(endIndex), '23*');
        for currIndex = endIndex - 1 : -1 : startIndex
            %  --O- 1 -
            %    |    |
            %    2    |
            %    |    |
            %  --O- 2 -
            temp = contract(psi(currIndex), 3, ERight, 1);
            ERight = contract(psi(currIndex), '23', temp, '23*');
        end
    else
        ERight = contract(psi(endIndex), '2', psi(endIndex), '2*', [1 3 2 4]);
        for currIndex = endIndex - 1 : -1 : startIndex
            %  --O-- 1 --|--
            %    |       |  
            %    2       | 
            %    |       |
            %  --O-- 2 --|--
            temp = contract(psi(currIndex), 3, ERight, 1);
            ERight = contract(temp, '23', psi(currIndex), '23*', [1 4 2 3]);
        end
    end
    E = ERight;
end