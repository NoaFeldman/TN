omegaOverTK = 1e-3;
% [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff1024, Gamma, omegaOverTK, 75, '', 1024);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = 1e-3; %TK * omegaOverTK * 1e-1;
for k = 2:3
    tic;
    [rAlphas, lAlphas, lambdaAlphas, rhoSS, Liou, rhoT] = diagLiou(NRG(k).HK, Sminus(k), SplusSminus(k), Inrg.EScale(k), gamma, T);
    omegas = Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 10;
    getOmegaAlphas = @(Sminus, rhoSS, rAlphas, lAlphas) getOmegaAlphasQSpace(Sminus, rhoSS, rAlphas, lAlphas);
    AsArnoldi = getSpectralFunction(omegas, Sminus(k), rhoSS, rAlphas, lAlphas, lambdaAlphas, getOmegaAlphas);
    % real(As)? abs(As)?
    [omegasSmooth, ASmooth] = getSmoothSpec(omegas, abs(AsArnoldi), 'emin', min(abs(omegas)), 'eps', T);
    % Normalize based on N from Scarlatllea's paper defined between eqs.
    % (5-6)
    ScarlatellaN = length(blkdiag(rhoT))/length(rAlphas);
    loglog(omegasSmooth, ASmooth/ScarlatellaN);
    hold on

    getOmegaAlphasM = @(Sminus, rhoSS, rAlphas, lAlphas) getOmegaAlphasMat(Sminus, rhoSS, rAlphas, lAlphas);
    toVec = @(QRho) reshape(blkdiag(QRho), [length(blkdiag(QRho))^2 1]);
    backToMat = @(vec) reshape(vec, [sqrt(length(vec)) sqrt(length(vec))]);

    rhoT = rhoT + 1e-99 * Liou.H;
    Liou.L = Liou.L + 1e-99 * Liou.H;
    Liou.LdagL = Liou.LdagL + 1e-99 * Liou.H;
    Sminus(k) = Sminus(k) + 1e-99 * Liou.H;
    rhoSS = rhoSS + 1e-50*Liou.H;
    id = getIdentity(Liou.H);
    applyLMatRight = @(Liou, r) blkdiag(Liou.H) * r - r * blkdiag(Liou.H) + ...
       blkdiag(Liou.L) * r * blkdiag(Liou.L') ...
       - 0.5 .* (blkdiag(Liou.LdagL) * r + r * blkdiag(Liou.LdagL));
    applyLMatLeft= @(Liou, l) blkdiag(Liou.H) * l - l * blkdiag(Liou.H) + ...
       blkdiag(Liou.L') * l * blkdiag(Liou.L) ...
       - 0.5 .* (blkdiag(Liou.LdagL) * l + l * blkdiag(Liou.LdagL));
    AsFull = zeros(length(omegas), 1);
    for b = 1:length(Liou.H.data)
        proj = cutQSpaceRows(id, b);
        LiouB = Liou;
        LiouB.H = LiouB.H*proj;
        LiouB.L = LiouB.L*proj;
        LiouB.LdagL = LiouB.LdagL*proj;
        rhoTB = rhoT * proj;
        rhoSSB = rhoSS * proj;
        SminusB = Sminus(k) * proj;
        [LiouMatB, rhoVecB] = vecing(LiouB, rhoTB);
        [evecsRB, evalsB, evecsLB] = eig(LiouMatB);
        evalsB = diag(evalsB);
        AsFull = AsFull + getSpectralFunction(omegas, blkdiag(SminusB), ...
            toVec(rhoSSB), evecsRB, evecsLB, evalsB, getOmegaAlphasM);
    end
    [omegasSmoothFull, ASmoothFull] = getSmoothSpec(omegas, abs(AsFull), ...
                                    'emin', min(abs(omegas)), 'eps', T);
    n = sqrt(length(blkdiag(rhoT)));
    ScarlatellaN = 1/n; % n / n^2
    loglog(omegasSmoothFull, ASmoothFull/ScarlatellaN);       
      %    
%     [LiouMatSub, basisSub, evecsRSub, evecsLSub, evalsSub] = ...
%         getSubBasis(rhoVec, LiouMat);
%     AsSub = getSpectralFunction(omegas, blkdiag(Sminus(k) + 1e-99 * NRG(k).HK), ...
%         toVec(rhoSS + 1e-50*rhoT), evecsRSub, evecsLSub, evalsSub, getOmegaAlphasM);
%     [omegasSmoothSub, ASmoothSub] = getSmoothSpec(omegas, abs(AsSub), ...
%                              'emin', min(abs(omegas)), 'eps', T);
%    loglog(omegasSmoothSub, ASmoothSub/ScarlatellaN);
%    
%    line([Inrg.EScale(k) Inrg.EScale(k)], [1e-9 1e-6], 'color', 'k');
%    legend({'Arnoldi', 'Exact - full basis', 'Exact - $\rho_T$ sub-basis', '$E_\mathrm{scale}$'}, 'interpreter', 'latex', 'fontsize', 20);
    
   toc;
   disp(strcat('k = ', int2str(k)));
end
