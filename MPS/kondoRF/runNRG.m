function [NRG, Inrg, AV, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, epsH, omegaL, Gamma, OmegaOverTK, N, fout)
    a = pi*(-U/(8*Gamma) + Gamma/(2*U));
    x = (epsE+U/2)*sqrt(pi/(2*U*Gamma));
    TK=min(1,sqrt(U*Gamma/2))*exp(a+x^2);
    Omega = OmegaOverTK * TK;
    % Get A0 and H0 for the Kondo problem with the valence and conduction
    % band as described in Sbierski et al
    [AV, A0, H0] = getH0(epsE, Ueh, U, epsH, omegaL, Omega);
    
    [F, Z, ~, ~] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    hackF = F;
    for i = 1:length(hackF)
        hackF(i).Q = hackF(i).Q(1:2);
        hackF(i).info.itags = hackF(i).info.itags(1:2);
    end
    
    % Calculate t_N values based on Bulla's paper
    Lambda = 2.7;
    ff = getNRGcoupling(Gamma, Lambda, N);
    
    [NRG, Inrg] = NRGWilsonQS(H0, A0, Lambda, ff, hackF, Z, 'fout', fout);
    NRG = makeNRGQSpace(NRG);
    EE=Inrg.EE;
    E0=Inrg.E0;
    
    AK = contract(AV, 2, A0, 1, [1 3 2 4]);
    id = getIdentity(AK, 3, AK, 4, 'sd');
    AK = contract(AK, '34', id, '12*');
    save(strcat(fout, '_00.mat'), 'AK', '-append');
    base = contract(id, 3, id, '3*');
    sigmaMinOp = contract(contract(base, 1, F(2)', 2, [4 1 2 3 5]), '25', F(2), '23', [1 4 2 3]);
    sigmaMinOp = contract(id, '12*', contract(sigmaMinOp, '34', id, '12'), '12');
end

function [NRGD] = makeNRGQSpace(NRGD)
    % Cast all NRGWilsonQS ouput operators to QSpace class
    for k = 1:length(NRGD)
        NRGD(k).AK = QSpace(NRGD(k).AK);
        NRGD(k).AT = QSpace(NRGD(k).AT);
        NRGD(k).HK = QSpace(NRGD(k).HK);
        NRGD(k).HT = QSpace(NRGD(k).HT);
    end
end