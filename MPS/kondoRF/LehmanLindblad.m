omegaOverTK = 1e-3;
[NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff1024, Gamma, omegaOverTK, 75, '', 1024);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = TK * omegaOverTK * 1e-1;
for k = 2:length(NRG)-1
   tic;
   [rAlphas, lAlphas, lambdaAplhas] = diagLiou(NRG(k).HK, Sminus(k), SplusSminus(k), Inrg.EScale(k), gamma, T);
   toc;
   disp(strcat('k = ', int2str(k)));
end
