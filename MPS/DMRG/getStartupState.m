function psi = getStartupState(N, m)
    % Create a product state (up-down-up-down-up...) for n sites with spin m.
    % Assuming n even.
    % Assuming m >= 0.
    % Bring the state to a left canonical form (except for the rightmost site).
    currSpin = 0;
    flipRate = N+1;
    if (m ~=0)
        flipRate = floor((N / 2) / m);
    end
    for i = 1:N
        psi(i) = QSpace;
        if (mod(i, 2) ~= 0 | mod(i / 2, flipRate) == 0)
            psi(i).Q = {[currSpin], [-1], [currSpin - 1]};
            currSpin = currSpin - 1;
        else
            psi(i).Q = {[currSpin], [+1], [currSpin  + 1]};
            currSpin = currSpin + 1;
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
    psi(N).info.itags{3} = strcat(int2str(N), 'a');
end
