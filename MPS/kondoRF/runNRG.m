function [NRGD, NRGI, AV] = runNRG(epsE, Ueh, U, epsH, omegaL, Gamma, N)
    Omega = 0;
    % Get A0 and H0 for the Kondo problem with the valence and conduction
    % band as described in Sbierski et al
    [AV, A0, H0] = getHO(epsE, Ueh, U, epsH, omegaL, Omega);
    
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    hackF = F;
    for i = 1:length(hackF)
        hackF(i).Q = hackF(i).Q(1:2);
        hackF(i).info.itags = hackF(i).info.itags(1:2);
    end
    
    % Calculate t_N values based on Bulla's paper
    Lambda = 2.7;
    tn = @(n) (1 + Lambda.^(-1)) .* (1 - Lambda.^(-n-1)) ./ ...
        (2 .* sqrt((1 - Lambda.^(-2.*n - 1)).*(1 - Lambda.^(-2.*n - 3))))...
        .*Lambda.^(-n/2);
    hopImpC0 = 0.5;
    ff = [hopImpC0 tn(0:N)];
    
    % Etrunc parameters identical to rnrg.m (I get the same problem if I
    % just don't add these flags)
    ETRUNC = [400 400 400 400 400 8];
    Etrunc = 5;
    [NRGD, NRGI] = NRGWilsonQS(H0, A0, Lambda, ff, hackF, Z, ...
        'ETRUNC', ETRUNC, 'Etrunc', Etrunc);
    NRGD = makeNRGQSpace(NRGD);
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