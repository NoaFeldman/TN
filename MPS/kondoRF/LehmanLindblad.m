omegaOverTK = 1e-3;
% [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff1024, Gamma, omegaOverTK, 75, '', 1024);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = 1e-3; %TK * omegaOverTK * 1e-1;
for k = 2
   tic;
   [rAlphas, lAlphas, lambdaAlphas, rhoSS, Liou, rhoT] = diagLiou(NRG(k).HK, Sminus(k), SplusSminus(k), Inrg.EScale(k), gamma, T);
   omegas = Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 10;
   getOmegaAlphas = @(Sminus, rhoSS, rAlphas, lAlphas) getOmegaAlphasQSpace(Sminus, rhoSS, rAlphas, lAlphas);
   AsArnoldi = getSpectralFunction(omegas, Sminus(k), rhoSS, rAlphas, lAlphas, lambdaAlphas, getOmegaAlphas);
   % real(As)? abs(As)?
   [omegasSmooth, ASmooth] = getSmoothSpec(omegas, real(AsArnoldi), 'emin', min(abs(omegas)), 'eps', T);
   % Normalize based on N from Scarlatllea's paper defined between eqs.
   % (5-6)
   ScarlatellaN = length(blkdiag(rhoT))/length(rAlphas);
   loglog(omegasSmooth, ASmooth/ScarlatellaN);
   hold on
   
   getOmegaAlphasM = @(Sminus, rhoSS, rAlphas, lAlphas) getOmegaAlphasMat(Sminus, rhoSS, rAlphas, lAlphas);
   toVec = @(QRho) reshape(blkdiag(QRho), [length(blkdiag(QRho))^2 1]);
   backToMat = @(vec) reshape(vec, [sqrt(length(vec)) sqrt(length(vec))]);
   [LiouMat, rhoVec, ~] = vecing(Liou, rhoT);
   n = sqrt(length(rhoVec));
   [evecsR, evals, evecsL] = eig(LiouMat);
   evals = diag(evals);
   Liou.L = Liou.L + 1e-99 * Liou.H;
   Liou.LdagL = Liou.LdagL + 1e-99 * Liou.H;
   applyLMatRight = @(Liou, r) blkdiag(Liou.H) * r - r * blkdiag(Liou.H) + ...
       blkdiag(Liou.L) * r * blkdiag(Liou.L') ...
       - 0.5 .* (blkdiag(Liou.LdagL) * r + r * blkdiag(Liou.LdagL));
   applyLMatLeft= @(Liou, l) blkdiag(Liou.H) * l - l * blkdiag(Liou.H) + ...
       blkdiag(Liou.L') * l * blkdiag(Liou.L) ...
       - 0.5 .* (blkdiag(Liou.LdagL) * l + l * blkdiag(Liou.LdagL));
   AsFull = getSpectralFunction(omegas, blkdiag(Sminus(k) + 1e-99 * NRG(k).HK), ...
       toVec(rhoSS + 1e-50*rhoT), evecsR, evecsL, evals, getOmegaAlphasM);
   [omegasSmoothFull, ASmoothFull] = getSmoothSpec(omegas, real(AsFull), ...
                                    'emin', min(abs(omegas)), 'eps', T);
   ScarlatellaN = 1/n; % n / n^2
   loglog(omegasSmoothFull, ASmoothFull/ScarlatellaN);
   
    [LiouMatSub, basisSub, evecsRSub, evecsLSub, evalsSub] = ...
        getDenseBasis(rhoVec, LiouMat);
    AsSub = getSpectralFunction(omegas, blkdiag(Sminus(k) + 1e-99 * NRG(k).HK), ...
        toVec(rhoSS + 1e-50*rhoT), evecsRSub, evecsLSub, evalsSub, getOmegaAlphasM);
    [omegasSmoothSub, ASmoothSub] = getSmoothSpec(omegas, abs(AsSub), ...
                             'emin', min(abs(omegas)), 'eps', T);
   loglog(omegasSmoothSub, ASmoothSub/ScarlatellaN);
    
   toc;
   disp(strcat('k = ', int2str(k)));
end
