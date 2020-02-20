function [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff, Gamma, OmegaOverTK, N, fout, Nkeep)
    % First option for TK is from the Bethe ansatz, second option is from
    % Sbierski et al's paper, they are in the same order of magnitude for
    % the parameters Sbierski et al (and so we) use.
    a = pi*(-U/(8*Gamma) + Gamma/(2*U));
    x = (epsE+U/2)*sqrt(pi/(2*U*Gamma));
    TK = min(1,sqrt(U*Gamma/2))*exp(a+x^2);
    TK = 1.5e-4;
    Omega = OmegaOverTK * TK;
    % Get A0 and H0 for the Kondo problem with the valence and conduction
    % band as described in Sbierski et al
    [AV, A0, H0, rabiOp] = getH0(epsE, Ueh, U, omegaDiff, Omega);
    
    [F, Z, ~, ~] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    hackF = F;
    for i = 1:length(hackF)
        hackF(i).Q = hackF(i).Q(1:2);
        hackF(i).info.itags = hackF(i).info.itags(1:2);
    end
    
    % Calculate t_N values based on Bulla's paper
    Lambda = 2.7;
    ff = getNRGcoupling(Gamma, Lambda, N);
    
    [NRG, Inrg] = NRGWilsonQS(H0, A0, Lambda, ff, hackF, Z, 'Nkeep', Nkeep, 'fout', fout);
    % Load NRG
    for i = 1:N
        NRG(i) = load(sprintf('%s_%02g.mat', fout, i-1));
    end
    % Cast to QSpace object
    NRG = makeNRGQSpace(NRG);
    EE=Inrg.EE;
    E0=Inrg.E0;
    AC = A0;
    AC.info.itags = {'L00', 'K00*', 's00'};
    
    % Include the valence level in the first NRG shell
    AK = contract(AV, 2, A0, 1, [1 3 2 4]);
    id = getIdentity(AK, 3, AK, 4, 's00');
    AK = contract(AK, '34', id, '12*');
    AK.info.itags{2} = 'K00*';
    save(strcat(fout, '_00.mat'), 'AK', '-append');
    NRG(1).AK = AK;
    sigmaMinOp = contract(id, '12*', contract(rabiOp, '34', id, '12'), '12');
    sigmaMinOp.info.itags = {'', '*'};
end

function [NRGD] = makeNRGQSpace(NRGD)
    % Cast all NRGWilsonQS ouput operators to QSpace class
    for k = 1:length(NRGD)
        NRGD(k).AK = QSpace(NRGD(k).AK);
        NRGD(k).AT = QSpace(NRGD(k).AT);
        NRGD(k).HK = diag(QSpace(NRGD(k).HK));
        NRGD(k).HT = diag(QSpace(NRGD(k).HT));
    end
end