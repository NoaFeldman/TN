function omegaAlphas = getOmegaAlphasQV(rAlphas, lDagAlphas, rhoSS, Sminus, Splus)
    % eq. (18) in https://arxiv.org/pdf/1811.03518.pdf.
    trRSminus = contract(rAlphas, '23', Sminus, '21');
    trRSminus = trRSminus.data{1};
    
    rhoSSSplus = rhoSS * Splus;
    SplusRhoSS = Splus * rhoSS;
    trLdagRhoSSSplus = contract(lDagAlphas, '23', rhoSSSplus, '21');
    trLdagRhoSSSplus = trLdagRhoSSSplus.data{1};
    trLdagSplusRhoSS = contract(lDagAlphas, '23', SplusRhoSS, '21');
    trLdagSplusRhoSS = trLdagSplusRhoSS.data{1};
    omegaAlphas = trRSminus .* (trLdagSplusRhoSS - trLdagRhoSSSplus);
end
