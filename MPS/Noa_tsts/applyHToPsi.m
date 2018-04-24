function Hpsi = applyHToPsi(H, psi)
    % TODO
    N = length(psi);
    
    % Apply H.single to all sites and sum all arrays to get
    % single x 1 x 1 ... + 1 x single x 1 ... + 1 x 1 x single ... + ...
    % save all operators for two sites separately, later when contracting
    % the sites to one another we will contract these as well.
    % TODO this win't work
    HPsi = QSpace(N);
    HR2L = QSpace(N-1);
    for i = 1 : N
        HSingle = psi;
        HSingle(i) = contract(H.single(i), 1, psi(i), 2, [2 1 3]);
        HPsi = HPsi + HSingle;
        if (i ~= N)
            HR2L(i) = contract( ...
                contract(psi(i), 3, psi(i+1), 1), '23', ...
                contract(H.l2r(i), 3, H.r2l(i+1), 3), '13', [1 3 4 2]);
        end
    end
    % Contract all sites
end
        
    