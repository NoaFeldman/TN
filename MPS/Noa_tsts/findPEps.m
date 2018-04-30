function res = findPEps(t, p)
    K = 1;
    n = 1;
    NA = -2:0;
    cn = 1;
    eps = 1:1e-1:20;
    tMat = repmat(t.', 1, length(eps));
    epsMat = repmat(eps, length(t), 1);
    x = (epsMat.^2 + tMat.^2) ./ (epsMat / 2);
    chi = zeros(1, length(eps));
    for na = 1:length(NA)
        pTheo = zeros(length(t), length(eps));
%         dalpha = 1e-4;
%         for alpha = -pi:dalpha:pi
%             pTheo = pTheo + cn .* x.^(-2 * (alpha/(2*pi))^2) * dalpha / (2*pi) ...
%                 .* exp(1i * alpha * NA(na));
%         end
        pTheo = cn .* sqrt(pi * n ./ (2 * K .*  log(x))) ...
        .* exp(-pi^2 * n * NA(na)^2 ./ (2 * K .* log(x)));
        res.(strcat('ptheo', int2str(na))) = pTheo;
        pMat = repmat(p(na).', 1, length(eps));
        chi = chi + (sum((pMat - pTheo).^2, 1));
        res.(strcat('chi', int2str(na))) = chi;
     end
end