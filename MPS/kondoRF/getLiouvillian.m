function Liou  = getLiouvillian(H, g, sMinus, sPlusSMinus)
    Liou = struct();
    Liou.H = H;
    Liou.L = sqrt(g) * sMinus;
    Liou.LdagL = g*sPlusSMinus;
end