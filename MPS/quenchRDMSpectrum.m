function quenchRDMSpectrum(L, h, JPM, JZ, m, dt, tStep, tFirstStep, tStepsNum, dirName)
    % Saves the eigenvalues of the RDM (of half the lattice) for two ground
    % state L/2 lattices suddenly coupled.
    % We first calculate g.s of an L/2 lattice.
    % Then we couple two chains in g.s to an L sized lattice and propegate
    % this lattice in time (With trotter time step = dt).
    % Every t = tStep * dt we save the spectrum of the RDM.
    % h, JPM, JZ - Heisenberg hamiltonian parameters.
    % m = initial spin of state (assumed spin of g.s) for the L / 2 chain.
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    tic;
    gsopts = {'Nkeep', 1024, 'stol', 1e-16};
    topts = {'Nkeep', 4};
%     dirName = strcat('quenchSpecL', int2str(L), 'JPM', num2str(abs(JPM)), ...
%                      'JZ', num2str(abs(JZ)), 'h', num2str(abs(h)), '_200');
    if (tFirstStep == 0)
        [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m, gsopts);
        disp('Found ground state for L/2');
        toc;
        [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
        psi = coupleStates(gs, gs);
        mkdir(dirName);
    else
        f = load(strcat(dirName, '/psiAtStep', int2str(tFirstStep - 1), '.mat'));
        psi = f.psi;
        [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
    end
    truncErr = zeros(1, tStep*(tStepsNum - tFirstStep + 1));
    trotterGates = getTrotterGates(H, dt, 0);
    for step = tFirstStep : tStepsNum
        saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), psi, length(psi)/4);
        if (step < tStepsNum)
            for t = 1 : tStep
                ind = (step - tFirstStep)*tStep + t;
                [psi, truncErr(ind)] = trotterSweep(trotterGates, psi, topts);
            end
        end
        if (mod(step, 200) == 0)
            save(strcat(dirName, '/psiAtStep', int2str(step), '.mat'), 'psi');
        end
    end
    save(strcat(dirName, '/truncErr.mat'), 'truncErr');
end