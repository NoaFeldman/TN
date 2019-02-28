function psi = getStartupState(N, m, bc)
    % Create a product state (up-down-up-down-up...) for n sites with spin m.
    % Assuming n even.
    % Assuming m >= 0.
    % Bring the state to a left canonical form (except for the rightmost site).
    if nargin == 1
        m = 0;
        bc = 'open'
    end
    if strcmp(bc, 'open')
%         currSpin = 0;
%         flipRate = N+1;
%         if (m ~=0)
%             flipRate = floor((N / 2) / m);
%         end
%         for i = 1:N
%             psi(i) = QSpace;
%             if (mod(i, 2) ~= 0 | mod(i / 2, flipRate) == 0)
%                 psi(i).Q = {[currSpin], [-1], [currSpin - 1]};
%                 currSpin = currSpin - 1;
%             else
%                 psi(i).Q = {[currSpin], [1], [currSpin + 1]};
%                 currSpin = currSpin + 1;
%             end
%             psi(i).data = {1};
%             tag1 = strcat(int2str(i - 1), 'a');
%             tag2 = strcat(int2str(i), 's');
%             tag3 = strcat(int2str(i), 'a', '*');
%             psi(i).info.itags = {tag1, tag2, tag3};
%             psi(i).info.qtype = '';
%             psi(i).info.otype = '';
%         end
%         % Get the rightmost site to mixed canonical form.
%         psi(N).Q{3} = -1 * psi(N).Q{3};
%         psi(N).info.itags{3} = strcat(int2str(N), 'a');
%         
        fullState = QSpace();
        fullState.info.qtype = '';
        fullState.info.otype = '';
        fullState.Q{1} = 0;
        fullState.info.itags{1} = '0a';
        for i = 2:N+1
            if (mod(i, 2) == 0)
                fullState.Q{i} = 1;
            else
                fullState.Q{i} = -1;
            end
            fullState.info.itags{i} = strcat(int2str(i-1), 's');
        end
        fullState.Q{N+2} = 0;
        fullState.info.itags{N+2} = strcat(int2str(N), 'a');
        fullState.data{1} = 1;
        
        fullStateFlipped = fullState;
        for i = 2:N+1
            fullStateFlipped.Q{i} = -fullState.Q{i};
        end
        fullState = fullState - fullStateFlipped;
        psi = QSpace(N);
        for i = 1:N-1
            [psi(i), fullState, ~] = orthoQS(fullState, [1 2], '>>', ...
                            'Nkeep', 1024, 'itag', strcat(int2str(i), 'a'));
        end
        psi(N) = fullState;
        norm = getOverlap(psi, psi);
        for i = 1:length(psi(1).data)
            psi(1).data{i} = psi(1).data{i}/sqrt(norm);
        end
        
        if m == -1
            % for now only working for |m|, if we go to the asymmetric case
            % fix.
            [S,IS]=getLocalSpace('Spin',0.5,'-A');
            S(2).Q{3} = -S(2).Q{3};
            S(2).info.itags{3} = '';
            S(2).info.itags = {strcat(int2str(N), 's'), strcat(int2str(N), 's*'), 'temp'};
            temp = contract(psi(N), 2, S(2), 2);
            id = getIdentity(temp, 2, temp, 4, psi(N).info.itags{3});
            psi(N) = contract(temp, '24', id, '12*')
        end
        if m ~= 0 && m ~= -1
            disp('------------------------Unsupported m! Can be fixed in getStartupState.------------------');
        end
    else
        Q = 0;
        for i = 1:N
            psi(i)= QSpace();
            if mod(i, 2) == 0
                psi(i).Q = {Q; 1; -1; Q};
            else
                psi(i).Q = {Q; -1; 1; Q};
            end
            psi(i).data = {1};
            psi(i).info.qtype = '';
            psi(i).info.otype = '';
            psi(i).info.itags = ...
                {strcat(int2str(i-1), 'a'), strcat(int2str(i), 'sUp'), ...
                strcat(int2str(i), 'sDown'), strcat(int2str(i),'a*')};
        end
        psi(N).info.itags{4} = strcat(int2str(N), 'a');
        psi(N).Q{4} = -psi(N).Q{4};
    end
end
