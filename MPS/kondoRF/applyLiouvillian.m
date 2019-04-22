function LRho = applyLiouvillian(Rho, Liou, dag)
    if ~exist('dag','var'), dag=0; end
    LRho = (-1)^dag * -1j * (Liou.H * Rho - Rho * Liou.H);
    for Li = 1:length(Liou.L)
        if dag, L = Liou.L(Li)'; else, L = Liou.L(Li); end
        LdagL = Liou.LdagL(Li);
        LRho = LRho + L*Rho*L' - 0.5*LdagL*Rho -0.5*Rho*LdagL;
    end
end