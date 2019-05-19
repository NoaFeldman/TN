load('/home/noa/TN/MPS/kondoRF/params.mat')
omegaOverTK = 1e-3;
fout = '/home/noa/TN/MPS/kondoRF/nrgRes';
N = 75;
Nkeep = 1024;
omegaDiff = omegaDiff1024;
% [NRG, Inrg, AV, AC, EE, E0, TK, sigmaMinOp] = runNRG(epsE, Ueh, U, omegaDiff, Gamma, omegaOverTK, N, fout, Nkeep);
T = 1e-11;
Sminus = getOp(NRG, sigmaMinOp);
Splus = getOp(NRG, sigmaMinOp');
SplusSminus = getOp(NRG, sigmaMinOp'*sigmaMinOp);
gamma = 0; %1e-3; %TK * omegaOverTK * 1e-1;
k = 10;
H = NRG(k).HK;
Sminus(k) = Sminus(k) + 1e-99 * H;
Splus(k) = Splus(k) + 1e-99 * H;
Liou = getLiouvillian(H, gamma, Sminus(k), SplusSminus(k));
rhoT = getThermalState(H, 1);
calcOverlap = @(rho1, rho2) trace(rho1'*rho2);
rhoT = rhoT / sqrt(calcOverlap(rhoT, rhoT));
% [om, a0, Idma] = fdmNRG_QS(...
%     [fout '/NRG'], sigmaMinOp, sigmaMinOp, QSpace(), 'emin', 1e-10, 'T', T);
% [omegasSmooth_FDM, ASmooth_FDM] = getSmoothSpec(om, abs(a0), 'emin', min(abs(om)), 'eps', T);
loglog(omegasSmooth_FDM, ASmooth_FDM);
hold on

omegas = Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 1e-2 : Inrg.EScale(k) * 10;
As = zeros(1, length(omegas));
for b = 2:length(Liou.H.data)
    LiouB = Liou;
    LiouB.H = cutQSpaceRows(LiouB.H, b);
    LiouB.L = cutQSpaceRows(LiouB.L, b);
    LiouB.LdagL = cutQSpaceRows(LiouB.LdagL, b);
    rhoTB = cutQSpaceRows(rhoT, b);
    SminusB = cutQSpaceRows(Sminus(k), b);
    SplusB = cutQSpaceRows(Splus(k), b);
    
    
    [vecedLiou, id, idIn] = qspaceVecing(LiouB, rhoTB);
    [rAlphas, lDagAlphas, lambdaAlphas, rhoSS] = diagVecedLiou(vecedLiou, id);
    omegaAlphas = getOmegaAlphasQV(rAlphas, lDagAlphas, rhoSS, SminusB, SplusB);
%     As = As + getSpectralFunctionQV(omegas, omegaAlphas, lambdaAlphas);
    As = getSpectralFunctionQV(omegas, omegaAlphas, lambdaAlphas);
    [omegasSmooth, ASmooth] = getSmoothSpec(omegas, abs(As)', 'emin', min(abs(omegas)), 'eps', T);
    loglog(omegasSmooth, ASmooth);
end

% Normalize!!!
