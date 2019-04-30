omegaOverTK = 1e-3;
[NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff1024, Gamma, omegaOverTK, 75, '', 1024);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
Splus = getOp(NRG, sigmaMinOp');
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = 1e-3; %TK * omegaOverTK * 1e-1;
k = 2;
H = NRG(k).HK;
Liou = getLiouvillian(H, gamma, Sminus(k), SplusSminus(k));
rhoT = getThermalState(H, 1);
calcOverlap = @(rho1, rho2) trace(rho1'*rho2);
rhoT = rhoT / sqrt(calcOverlap(rhoT, rhoT));

Liou2 = Liou;
Liou2.H = cutQSpaceRows(Liou2.H, 2);
Liou2.L = cutQSpaceRows(Liou2.L, 2);
Liou2.LdagL = cutQSpaceRows(Liou2.LdagL, 2);
rhoT2 = cutQSpaceRows(rhoT, 2);
[LiouMat, rhoVec] = vecing(Liou2, rhoT2);
[evecsR, evals, evecsL] = eig(LiouMat);
evals = diag(evals);
backToMat = @(vec) reshape(vec, [sqrt(length(vec)) sqrt(length(vec))]);
Sminus2 = cutQSpaceRows(Sminus(k) + 1e-99 * H, 2);
Splus2 = cutQSpaceRows(Splus2(k) + 1e-99 * H, 2);

for alpha = length(evals)
    trRs(alpha) = trace(backToMat(evecsR(:, i)) * Sminus2.data{1});
%     trLs(alpha) = 

end