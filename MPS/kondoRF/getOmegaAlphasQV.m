function omegaAlphas = getOmegaAlphasQV(rAlphas, lDagAlphas, Sminus, Splus, rhoSS, id, idIn)
    % eq. (18) in https://arxiv.org/pdf/1811.03518.pdf.
    SminusAdjusted =  contract(Sminus, 1, idIn, '1*');
    trRSminus = contract(rAlphas, '23', SminusAdjusted, '12');
    trRSminus = trRSminus.data{1};
    
    rhoSSSplus = idIn' * rhoSS * Splus;
    SplusRhoSS = idIn' * Splus * rhoSS;
    trLdagRhoSSSplus = contract(lDagAlphas, '23', rhoSSSplus, '12');
    trLdagRhoSSSplus = trLdagRhoSSSplus.data{1};
    trLdagSplusRhoSS = contract(lDagAlphas, '23', SplusRhoSS, '12');
    trLdagSplusRhoSS = trLdagSplusRhoSS.data{1};
    omegaAlphas = trRSminus .* (trLdagSplusRhoSS - trLdagRhoSSSplus);
end
