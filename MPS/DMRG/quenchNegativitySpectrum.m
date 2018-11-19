function quenchNegativitySpectrum(L, l, h, JPM, JZ, m, dt, tStep, tFirstStep, tStepsNum, dirName, cpuNumber) 
    % Saves the eigenvalues of the RDM (of half the lattice) for two ground
    % state L/2 lattices suddenly coupled.
    % We first calculate g.s of an L/2 lattice.
    % Then we couple two chains in g.s to an L sized lattice and propegate
    % this lattice in time (With trotter time step = dt).
    % Every t = tStep * dt we save the spectrum of the RDM.
    % h, JPM, JZ - Heisenberg hamiltonian parameters.
    % m = initial spin of state (assumed spin of g.s) for the L / 2 chain.
    maxNumCompThreads(cpuNumber);
    path(path, ['/home/noa/MPS/MPSPACK_v3.0']);
    path(path, ['/home/noa/MPS/DMRG']);
    startup;
    tic;
    topts = {'Nkeep', 1024, 'stol', 1e-8};
%     dirName = strcat('quenchSpecL', int2str(L), 'JPM', num2str(abs(JPM)), ...
%                      'JZ', num2str(abs(JZ)), 'h', num2str(abs(h)), '_200');
    if (tFirstStep == 0)
        [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m);
        disp('Found ground state for L/2');
        toc;
        [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
        trotterGates = getTrotterGates(H, dt, 0);
        mkdir(dirName);
        save(strcat(dirName, '/trotterGates.mat'), 'trotterGates');
        psi = coupleStates(gs, mirrorState(gs));
        saveNegativitySpectrum(psi, l, 'symm', strcat(dirName, '/step0'));
    else
        f = load(strcat(dirName, '/psiAtStep', int2str(tFirstStep - 1), '.mat'));
        psi = f.psi;
        f = load(strcat(dirName, '/trotterGates.mat'));
        trotterGates = f.trotterGates;
%         [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
    end
    truncErr = zeros(1, tStep*(tStepsNum - tFirstStep + 1));
    
    t = tFirstStep;
    vF = 2;
    for step = tFirstStep : tStepsNum
        for tInd = 1 : tStep
            t = t + dt;
            ind = (step - tFirstStep)*tStep + tInd;
            [psi, truncErr(ind)] = trotterSweep(trotterGates, psi, topts);
        end
        saveNegativitySpectrum(psi, l, 'symm', strcat(dirName, '/step', int2str(step + 1)));
        if (mod(step, 1) == 0)
            save(strcat(dirName, '/psiAtStep', int2str(step), '.mat'), 'psi', 'truncErr');
        end
    end
    save(strcat(dirName, '/truncErr.mat'), 'truncErr');
end