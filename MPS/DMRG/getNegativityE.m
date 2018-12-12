function E = getNegativityE(psi, startIndex, endIndex, dir)
    % compute E matrices from figure 6 at https://arxiv.org/pdf/1605.00674.pdf
    % Working site is expected to be psi(endIndex).
    % startIndex and endIndex start from the right.
    if dir == '<<'
        ERight = contract(QSpace(psi(endIndex)), '2', QSpace(psi(endIndex)), '2*');
        for currIndex = endIndex - 1 : -1 : startIndex
            %  --O-- 1 --|--
            %    |       |  
            %    2       | 
            %    |       |
            %  --O-- 2 --|--
            temp = contract(QSpace(psi(currIndex)), 3, ERight, 1);
            ERight = contract(temp, '24', psi(currIndex), '23*', [1 2 4 3]);
            w = whos('ERight');
            disp(strcat('currIndex = ', int2str(currIndex), ', size = ', num2str(w.bytes * 1e-9), 'GB'));
        end
        E = ERight;
    else
        ELeft = contract(QSpace(psi(startIndex)), '2', QSpace(psi(startIndex)), '2*');
        for currIndex = startIndex + 1 : endIndex
            %  --|-- 1 --O--
            %    |       |  
            %    |       2 
            %    |       |
            %  --|-- 2---O--
            temp = contract(ELeft, 2, QSpace(psi(currIndex)), 1);
            ELeft = contract(temp, '34', psi(currIndex), '12*', [1 3 2 4]);
            w = whos('ELeft');
            disp(strcat('currIndex = ', int2str(currIndex), ', size = ', num2str(w.bytes * 1e-9), 'GB'));
        end
        E = ELeft;
    end
end