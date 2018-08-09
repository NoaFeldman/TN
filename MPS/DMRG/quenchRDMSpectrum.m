function quenchRDMSpectrum(L, h, JPM, JZ, m, dt, tStep, tFirstStep, tStepsNum, dirName, subsystemA)
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
    topts = {'Nkeep', 100};
%     dirName = strcat('quenchSpecL', int2str(L), 'JPM', num2str(abs(JPM)), ...
%                      'JZ', num2str(abs(JZ)), 'h', num2str(abs(h)), '_200');
    if (tFirstStep == 0)
        [gs, ~, ~, ~] = getGroundState(L / 2, h, JPM, JZ, m);
        disp('Found ground state for L/2');
        toc;
        [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
        psi = coupleStates(gs, mirrorState(gs));
        mkdir(dirName);
    else
        f = load(strcat(dirName, '/psiAtStep', int2str(tFirstStep - 1), '.mat'));
        psi = f.psi;
        [~, H, ~, ~] = myStartup(L, h, JPM, JZ, m);
    end
    truncErr = zeros(1, tStep*(tStepsNum - tFirstStep + 1));
    trotterGates = getTrotterGates(H, dt, 0);
    if nargin == 10
        subsystemA = 'half';
    end
    for step = tFirstStep : tStepsNum
        if strcmp(subsystemA, 'half')
            saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), psi, length(psi)/2);
        % For the lightcone states, we take the B system to be ~L/16, so
        % the plateau will be short.
        elseif strcmp(subsystemA, 'asymEdge')
            saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), psi, 7 * length(psi)/16);
        elseif strcmp(subsystemA, 'symMid')
            folded = foldChain(psi);
            saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), folded, length(folded)/8);
        elseif strcmp(subsystemA, 'asymMid')
            aux = getStartupState(7 * L / 16);
            % [aux   - B1   - A     - B2]
            % [7L/16 - L/16 - 7L/16 - 8L/16]
            % We fold in the middle of subsystem A:
            % [aux - B1 # A/2]
            % [   B2    # A/2]
            % We measure the entanglement in #: 8L/16
            coupled = coupleStates(aux, psi);
            folded = foldChain(coupled);
            saveRDMSpectrum(strcat(dirName, '/step', int2str(step), '.mat'), folded, L/2);
        end
        
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