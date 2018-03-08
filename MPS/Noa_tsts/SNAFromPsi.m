function SNAFromPsi(psi, spectrumFileName, figFileName)
    saveRDMSpectrum(spectrumFileName, psi);
    spectrum = load(spectrumFileName);
    sz = keys(spectrum);
    val = values(spectrum);
    for i = 1:length(sz)
        x(i) = str2num(sz{i});
        s(i) = -sum(val{i}.*log(val{i}));
    end
    STheo = sqrt(pi * log(length(psi)) / 2) / 3 .* exp(-1 * pi^2 .* (x/2).^2 ./ (2 * log(length(psi))));
    scatter(x(:), s(:), 'markerEdgeColor', 'r')
    hold on
    scatter(x(:), STheo(:), 'markerEdgeColor', 'g')
    legend({'DMRG', 'Theoretical'})
    xlabel('2$S^Z_A$', 'Interpreter', 'latex');
    ylabel('$S(N_A)$', 'Interpreter', 'latex');
    savefig(figFileName);
end