function omegaAlphas = getOmegaAlphasQV(rAlphas, lDagAlphas, rhoSS, Sminus, Splus)
    % eq. (18) in https://arxiv.org/pdf/1811.03518.pdf.
    
    %  s^- * r_alpha:
    %
    %    k1
    %    |
    %   |s|   k2
    %    |     |
    %    |__r__|
    %       | 
    %     alpha
    % And we trace this expression by contracting the k1, k2 legs. we are
    % left with a rank 1 tensor of all the alpha values.
    trRSminus = contract(rAlphas, '12', Sminus, '21');
    trRSminus = trRSminus.data{1};
    
    % Same trick for Tr(l^dagger_alpha * [s^+, rho_s])
    rhoSSSplus = rhoSS * Splus;
    SplusRhoSS = Splus * rhoSS;
    trLdagRhoSSSplus = contract(lDagAlphas, '12', rhoSSSplus, '21');
    trLdagRhoSSSplus = trLdagRhoSSSplus.data{1};
    trLdagSplusRhoSS = contract(lDagAlphas, '12', SplusRhoSS, '21');
    trLdagSplusRhoSS = trLdagSplusRhoSS.data{1};
    
    omegaAlphas = trRSminus .* (trLdagSplusRhoSS - trLdagRhoSSSplus);
end
