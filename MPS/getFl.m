function getFl(L, h, JPM, JZ, m, dt, tStep, tStepsNum)
    % Creates vectors of fl, the eigenvalues of <Ci^dagger * Cj> (based on
    % page 2 of https://arxiv.org/pdf/1711.09418.pdf).
    % We first calculate g.s of an L/2 lattice.
    % Then we couple two chains in g.s to an L sized lattice and propegate
    % this lattice in time (With trotter time step = dt).
    % Every t = tStep * dt we calculate <Ci^dagger * Cj> and diagonalize.
    % h, JPM, JZ - Heisenberg hamiltonian parameters.
    % m = initial spin of state (assumed spin of g.s) for the L / 2 chain.
    opts = {'Nkeep', 2048, 'stol', 1e-5};
    [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m, opts);
    [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
    psi = coupleStates(gs, gs);
    for step = 1: tStepsNum
        for t = 1 : tStep
           psi = trotterSweep(psi, dt, 0, H, opts);
        end
        M = getCiCjMatrix(psi, L / 2);
        [~, F] = eig(M);
        f = zeros(1, length(F));
        for i = 1 : length(F)
            f(i) = F(i, i);
        end
        res.(strcat('flT', num2str(step))) = f;
    end
    save(strcat('fL', int2str(L), 'JPM', abs(num2str(JPM)), ...
            'JZ', abs(num2str(JZ)), 'h', abs(num2str(h)), '.mat'), ...
        'res');
end
    