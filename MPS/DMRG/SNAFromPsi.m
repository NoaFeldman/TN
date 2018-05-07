function SNAFromPsi(psi, spectrumFileName, figFileName)
    saveRDMSpectrum(spectrumFileName, psi);
    s = load(spectrumFileName);
    sz = keys(s.spectrum);
    val = values(s.spectrum);
    x = zeros(1, length(sz));
    s = zeros(1, length(sz));
    for i = 1:length(sz)
        x(i) = str2num(sz{i});
%         s(i) = -sum(val{i}.*log2(val{i}));
        dn = 1e-2;
        s(i) = -(sum(val{i}.^(1+dn)) - sum(val{i}.^(1-dn))) / (2 * dn);
    end
    [x, x_order] = sort(x);
    s = s(x_order);
    STheo = sqrt(pi * log(length(psi)) / 2) / 3 .* exp(-1 * pi^2 .* (x/2).^2 ./ (2 * log(length(psi))));
    scatter(x(:), s(:), 'markerEdgeColor', 'r');
    hold on
    plot(x, STheo, 'color', 'r');
%     legend({'DMRG', 'Theoretical'});
    xlabel('2$S^Z_A$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    savefig(figFileName);
    hold off
end