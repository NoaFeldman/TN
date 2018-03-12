function quenchRDMSpectrum(L, h, JPM, JZ, m, dt, tStep, tStepsNum)
    % Saves the eigenvalues of the RDM (of half the lattice) for two ground
    % state L/2 lattices suddenly coupled.
    % We first calculate g.s of an L/2 lattice.
    % Then we couple two chains in g.s to an L sized lattice and propegate
    % this lattice in time (With trotter time step = dt).
    % Every t = tStep * dt we save the spectrum of the RDM.
    % h, JPM, JZ - Heisenberg hamiltonian parameters.
    % m = initial spin of state (assumed spin of g.s) for the L / 2 chain.
    tic;
    opts = {'Nkeep', 2048, 'stol', 1e-5};
    [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m, opts);
    disp('Found ground state for L/2');
    toc;
    [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
    psi = coupleStates(gs, gs);
    dirName = strcat('quenchSpecL', int2str(L), 'JPM', num2str(abs(JPM)), ...
                     'JZ', num2str(abs(JZ)), 'h', num2str(abs(h)));
    mkdir(dirName);
    truncErr = 0;
    for step = 0: tStepsNum
        saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), psi, truncErr);
        disp('saved RDM spectrum');
        toc;
        if (step < tStepsNum)
            for t = 1 : tStep
                [psi, truncErr] = trotterSweep(psi, dt, 0, H, opts);
                disp(strcat('t = ', int2str(t + (step-1)*tStep) , ' * ', num2str(dt)));
                disp(strcat('truncErr = ', num2str(truncErr)));
                toc;
            end
        end
    end
end