function Liou  = getLiouvillian(H, g, sMinus, sPlusSMinus)
    Liou = struct();
    Liou.H = H;
    Liou.L = sqrt(g) * sMinus + 1e-99 * H;
    Liou.LdagL = g*sPlusSMinus + 1e-99 * H;
end