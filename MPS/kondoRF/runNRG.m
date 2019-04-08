function [NRGD, NRGI, AV, EE, E0] = runNRG(epsE, Ueh, U, epsH, omegaL, Gamma, N, fout)
    a = pi*(-U/(8*Gamma) + Gamma/(2*U));
    x = (epsE+U/2)*sqrt(pi/(2*U*Gamma));
    TK=min(1,sqrt(U*Gamma/2))*exp(a+x^2);
    Omega = 2e-10;
    % Get A0 and H0 for the Kondo problem with the valence and conduction
    % band as described in Sbierski et al
    [AV, A0, H0] = getH0(epsE, Ueh, U, epsH, omegaL, Omega);
    
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    hackF = F;
    for i = 1:length(hackF)
        hackF(i).Q = hackF(i).Q(1:2);
        hackF(i).info.itags = hackF(i).info.itags(1:2);
    end
    
    % Calculate t_N values based on Bulla's paper
    Lambda = 2.7;
    ff = getNRGcoupling(Gamma, Lambda, N);
    
    [NRGD, NRGI] = NRGWilsonQS(H0, A0, Lambda, ff, hackF, Z); %, 'fout', fout);
    NRGD = makeNRGQSpace(NRGD);
    EE=NRGI.EE;
    E0=NRGI.E0;
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