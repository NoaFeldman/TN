function psi = getStartupState(N)
    % Create a product state (up-down-up-down-up...) for n sites.
    % Assuming n even.
    % Bring the state to a left canonical form (except for the rightmost site).
    for i = 1:N
        psi(i) = QSpace;
        if (mod(i, 2) == 0)
            psi(i).Q = {[-1], [-1], [0]};
        else
            psi(i).Q = {[0], [1], [-1]};
        end
        psi(i).data = {1};
        tag1 = strcat(int2str(i - 1), 'a');
        tag2 = strcat(int2str(i), 's');
        tag3 = strcat(int2str(i), 'a', '*');
        psi(i).info.itags = {tag1, tag2, tag3};
        psi(i).info.qtype = '';
        psi(i).info.otype = '';
    end
    % Get the rightmost site to mixed canonical form.
    psi(N).Q{3} = -1 * psi(N).Q{3};
    psi(N).info.itags{3} = strcat('a', int2str(N));
    
%     % chop the legs of the rightmost and leftmost sites
%     psi(1).Q = psi(1).Q(2:3);
%     psi(1).info.itags = psi(1).info.itags(2:3);
%     psi(N).Q = psi(N).Q(1:2);
%     psi(N).info.itags = psi(N).info.itags(1:2);
