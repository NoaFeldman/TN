function chargeConservedExact(paramFileName, fout, N, Nkeep, gamma, ...
    artificiallySmooth, compareToClosedSystem, outDirName, outFigureName, ...
    basis)
    % Notations and calculation are based on https://arxiv.org/pdf/1811.03518.pdf
    % artificiallySmooth - 0 or 1 variable that decides whether we add a
    % small real part to lambdaAlphas (eigenvalues of the Liouvillian).
    % compareToClosedSystem - Compare spectral functions to closed system
    % Lehman functions, for gamma = 0
    % basis - 'Full' for KK, KT, TK, TT. 'KK' for only KK

    mkdir outDirName
    load(paramFileName); % '/home/noa/TN/MPS/kondoRF/params.mat')
    omegaOverTK = 1e-3;
    omegaDiff = omegaDiff1024;
    [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff, Gamma, omegaOverTK, N, fout, Nkeep);
    T = Inrg.EScale(N);
    SminusKKArr = getOpKK(NRG, sigmaMinOp);
    SplusKKArr = getOpKK(NRG, sigmaMinOp');
    SplusSminusKKArr = getOpKK(NRG, sigmaMinOp'*sigmaMinOp);

    [om, a0, Idma] = fdmNRG_QS(...
        fout, sigmaMinOp, sigmaMinOp, QSpace(), 'emin', 1e-10, 'T', T);
    [omegasSmooth_FDM, ASmooth_FDM] = getSmoothSpec(om, abs(a0), 'emin', min(abs(om)), 'eps', T);

    omegas = 1.1.^(-240:25); % 1.1 .^ (~log(1.1, 1e-10):~log(1.1, 10))
    % TODO maybe adjust omegas for each shell separately?
    % NOTE: If artificiallySmooth is on, omega should have a much narrower
    % resolution.
    maxK = N;

    tic;
    for k = 1:maxK
        As = zeros(1, length(omegas));
        AsClosed = zeros(1, length(omegas));
        if strcmp(basis, 'Full')
            H = RabiKondo.getNRGfullBasis(diag(NRG(k).HK)', diag(NRG(k).HT)');
            H.data = cellfun(@diag, H.data, 'uniformOutput', 0);
            H.info.itags{1}(1) = 'F';
            AF = RabiKondo.getNRGfullBasis(NRG(k).AK, NRG(k).AT);

            Sminus = getOpInFullNRGBasis(SminusKKArr, k, NRG(k), H);
            Splus = getOpInFullNRGBasis(SplusKKArr, k, NRG(k), H);
            SplusSminus = getOpInFullNRGBasis(SplusSminusKKArr, k, NRG(k), H);
        elseif strcmp(basis, 'KK')
            if k == N
                H = NRG(k).HT;
            else
                H = NRG(k).HK;
            end
            Sminus = SminusKKArr(k);
            Splus = SplusKKArr(k);
            SplusSminus = SplusSminusKKArr(k);
        end

        Liou = getLiouvillian(H, gamma, Sminus, SplusSminus);
        rhoT = getThermalState(H, T);
        rhoT = rhoT + 1e-99 * H;
        rhoT = rhoT / trace(rhoT);

        for b = 11:length(Liou.H.data)
            LiouB = Liou;
            LiouB.H = cutQSpaceRows(LiouB.H, b);
            LiouB.L = cutQSpaceRows(LiouB.L, b);
            LiouB.LdagL = cutQSpaceRows(LiouB.LdagL, b);
            rhoTB = cutQSpaceRows(rhoT, b);
            weightB = trace(rhoTB);
            SminusB = cutQSpaceRows(Sminus, b);
            SplusB = cutQSpaceRows(Splus, b);

            [vecedLiou, id, idIn] = vecing(LiouB, LiouB.H);
            [rAlphas, lAlphas, lambdaAlphas, rhoSS] = diagVecedLiou(vecedLiou, id);
            if artificiallySmooth
                % Artificially smooth the spectrum by ensuring lambda_alpha has a
                % real value
                % NOTE: This option requires a much narrower resolution on
                % omega
                lambdaAlphas = lambdaAlphas - 1e-15 .* abs(imag(lambdaAlphas));
            end
            omegaAlphas = getOmegaAlphas(rAlphas, lAlphas, rhoSS, SminusB, SplusB);
            AsB =  getSpectralFunction(omegas, omegaAlphas, lambdaAlphas);
            As = As + weightB .* AsB;

            if compareToClosedSystem            
                AsClosedB =  zeros(1, length(omegas));
                % rhoT is diagonalized in the basis we work in
                % Eq. 16 in https://arxiv.org/pdf/1811.03518.pdf
                for n = 1:length(rhoTB.data{1})
                    for m = 1:length(rhoTB.data{1})
                        [~, ind] = min(abs(omegas - ...
                                   (LiouB.H.data{1}(m, m) - LiouB.H.data{1}(n, n))));
                        AsClosedB(ind) = AsClosedB(ind) + ...
                            (rhoTB.data{1}(n, n) - rhoTB.data{1}(m, m)) * ...
                                abs(SminusB.data{1}(n, m))^2;
                    end
                end
                AsClosed = AsClosed + AsClosedB;           
            end
        end
        disp(['finished shell ' int2str(k)]);
        toc;
        save([outDirName '/shell' int2str(k)], 'As', 'AsClosed', 'omegas');
    end
    
    fig = figure('visible','off');
    loglog(omegasSmooth_FDM/TK, ASmooth_FDM*TK);
    legendInfo = {'FDM'};
    hold on
    for k = 2:2:maxK-1
        a = load([outDirName '/shell' int2str(k)]);
        a2 = load([outDirName '/shell' int2str(k+1)]);
        [omegasSmooth, ASmooth] = ...
            getSmoothSpec(a.omegas, abs((a.As + a2.As) / 2)', ...
                          'emin', min(abs(a.omegas)), 'eps', T);
        loglog(omegasSmooth, ASmooth, 'color', [1/maxK*k 0 0.4]);
        legendInfo{end+1} = ['k' int2str(k)];
    end
    xlabel('$\omega$', 'interpreter', 'latex');
    ylabel('$A(\omega)$', 'interpreter', 'latex');
    legend(legendInfo);
    saveas(fig, [outDirName '/' outFigureName], 'fig')
    
    hold off
    if compareToClosedSystem
        fig = figure('visible','off');
        loglog(omegasSmooth_FDM/TK, ASmooth_FDM*TK);
        legendInfo = {'FDM'};
        hold on
        for k = 2:2:maxK-1
            a = load([outDirName '/shell' int2str(k)]);
            a2 = load([outDirName '/shell' int2str(k+1)]);
            [omegasSmooth, ASmooth] = ...
            getSmoothSpec(a.omegas, abs((a.AsClosed + a2.AsClosed) / 2)', ...
                          'emin', min(abs(a.omegas)), 'eps', T);
            loglog(omegasSmooth, ASmooth, 'color', [1/maxK*k 0 0.4]);
            legendInfo{end+1} = ['k' int2str(k)];
        end
        xlabel('$\omega$', 'interpreter', 'latex');
        ylabel('$A(\omega)$', 'interpreter', 'latex');
    end
    xlabel('$\omega$', 'interpreter', 'latex');
    ylabel('$A(\omega)$', 'interpreter', 'latex');
    legend(legendInfo);
    saveas(fig, [outDirName '/' outFigureName '_closed'], 'fig')
end
    
        
