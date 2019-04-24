omegaOverTK = 1e-3;
% [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff1024, Gamma, omegaOverTK, 75, '', 1024);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = 1e-3; %TK * omegaOverTK * 1e-1;
for k = 2:length(NRG)-1
   tic;
   [rAlphas, lAlphas, lambdaAlphas, rhoSS, Liou, rhoT] = diagLiou(NRG(k).HK, Sminus(k), SplusSminus(k), Inrg.EScale(k), gamma, T);
   omegas = Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 10;
   getOmegaAlphas = @(Sminus, rhoSS, rAlphas, lAlphas) getOmegaAlphasQSpace(Sminus, rhoSS, rAlphas, lAlphas);
   As = getSpectralFunction(omegas, Sminus(k), rhoSS, rAlphas, lAlphas, lambdaAlphas, getOmegaAlphas);
   [omegasSmooth,ASmooth] = getSmoothSpec(omegas, real(As), 'emin', min(abs(omegas)), 'eps', T);
   loglog(omegasSmooth, ASmooth);
   hold on
   toc;
   disp(strcat('k = ', int2str(k)));
end
