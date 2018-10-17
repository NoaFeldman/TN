function [I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, L2] = getRenyiNegNecessities(cicj, u1, v1, v2)
    u2 = v1+1;
    Gplus = getG(1, cicj, u1, v1, v2);
    Gminus = getG(-1, cicj, u1, v1, v2);
    I = eye (v2 - u1 + 1);
    Gp2 = Gplus^2;
    Gp3 = Gp2 * Gplus;
    GpGm = Gplus * Gminus;
    Gp2Gm = Gp2 * Gminus;
    L2 = v2 - u2 + 1;
end