function omegaAlphas = getOmegaAlphasQV(rAlphas, lDagAlphas, rhoSS, Sminus, Splus, id, idIn)
    % eq. (18) in https://arxiv.org/pdf/1811.03518.pdf.
    trRSminus = contract(rAlphas, '23', Sminus, '12');
    trRSminus = trRSminus.data{1};
    
    rhoSSSplus = idIn' * rhoSS * idIn' * Splus;
    SplusRhoSS = idIn' * Splus * rhoSS * idIn';
    trLdagRhoSSSplus = contract(lDagAlphas, '23', rhoSSSplus, '12');
    trLdagRhoSSSplus = trLdagRhoSSSplus.data{1};
    trLdagSplusRhoSS = contract(lDagAlphas, '23', SplusRhoSS, '12');
    trLdagSplusRhoSS = trLdagSplusRhoSS.data{1};
    omegaAlphas = trRSminus .* (trLdagSplusRhoSS - trLdagRhoSSSplus);
end
