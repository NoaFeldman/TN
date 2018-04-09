function quenchRDMSpectrum(L, h, JPM, JZ, m, dt, tStep, tFirstStep, tStepsNum)
    % Saves the eigenvalues of the RDM (of half the lattice) for two ground
    % state L/2 lattices suddenly coupled.
    % We first calculate g.s of an L/2 lattice.
    % Then we couple two chains in g.s to an L sized lattice and propegate
    % this lattice in time (With trotter time step = dt).
    % Every t = tStep * dt we save the spectrum of the RDM.
    % h, JPM, JZ - Heisenberg hamiltonian parameters.
    % m = initial spin of state (assumed spin of g.s) for the L / 2 chain.
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/Noa_tsts']);
    startup;
    tic;
    opts = {'Nkeep', 300};
    dirName = strcat('quenchSpecL', int2str(L), 'JPM', num2str(abs(JPM)), ...
                     'JZ', num2str(abs(JZ)), 'h', num2str(abs(h)));
    if (tFirstStep == 0)
        [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m, opts);
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
    for step = tFirstStep : tStepsNum
        saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), psi);
        disp('saved RDM spectrum');
        toc;
        if (step < tStepsNum)
            for t = 1 : tStep
                ind = (step - tFirstStep)*tStep + t;
                [psi, truncErr(ind)] = trotterSweep(psi, dt, 0, H, opts);
                disp(strcat('t = ', int2str(ind) , ' * ', num2str(dt)));
                disp(strcat('truncErr = ', num2str(truncErr(ind)))); % TODO relevant only for stol = 0
                toc;
            end
        end
        if (mod(step, 200) == 0)
            save(strcat(dirName, '/psiAtStep', int2str(step), '.mat'), 'psi');
        end
    end
    save(strcat(dirName, '/truncErr.mat'), 'truncErr');
end