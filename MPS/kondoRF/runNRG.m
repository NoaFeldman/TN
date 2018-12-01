function [NRGD, NRGI] = runNRG(epsE, Ueh, U, epsH, omegaL, Gamma, N)

    rho = 1/2;
    Tk = exp(-1/(Gamma / pi * (1 / abs(epsE) + 1/(epsH + U)) * 2));
    Omega = 0; % Tk * 1e-4;
    [AV, A0, H0] = getHO(epsE, Ueh, U, epsH, omegaL, Omega);
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    hackF = F;
    for i = 1:length(hackF)
    hackF(i).Q = hackF(i).Q(1:2);
    hackF(i).info.itags = hackF(i).info.itags(1:2);
    end
    Lambda = 2.7;
    tn = @(n) (1 + Lambda.^(-1)) .* (1 - Lambda.^(-n-1)) ./ (2 .* sqrt((1 - Lambda.^(-2.*n - 1)).*(1 - Lambda.^(-2.*n - 3)))).*Lambda.^(-n/2);
    hopImpC0 = sqrt(Gamma / (pi * rho));
    [NRGD, NRGI] = NRGWilsonQS(H0, A0, Lambda, [hopImpC0 tn(0:N)], hackF, Z);
    NRGD = makeNRGQSpace(NRGD);
end

function [NRG] = makeNRGQSpace(NRG)
    % Cast all NRGWilsonQS ouput operators to QSpace class
    for k = 1:length(NRG)
        NRG(k).AK = QSpace(NRG(k).AK);
        NRG(k).AT = QSpace(NRG(k).AT);
        NRG(k).HK = QSpace(NRG(k).HK);
        NRG(k).HT = QSpace(NRG(k).HT);
    end
end