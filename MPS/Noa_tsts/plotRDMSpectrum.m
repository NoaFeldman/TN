function plotRDMSpectrum(N, Delta)
    tic;
    % Recreate plot fFig.2 from https://arxiv.org/pdf/1407.3779.pdf
    [psi, H, HR, HL] = myStartup(N, 0, 1, Delta);
    % Find ground state
    ECurr = 0;
    opts = {'Nkeep', 100};
    for i=1:100
        EForm = ECurr;
        [HL, HR, psi, ~] = dmrgSweep(HL, HR, H, psi, '<<', opts);
        [HL, HR, psi, ECurr] = dmrgSweep(HL, HR, H, psi, '>>', opts);
        if (abs(ECurr - EForm)/abs(ECurr) < 1e-5)
            break;
        end
    end
    % Sweep to mid chain 
    i = length(psi);
    while(i > floor(length(psi)/2))
        [HL, HR, psi, E, i, M] = dmrgStep(HL, HR, H, psi, i, '<<', opts);
    end
    [~, ~, I] = orthoQS(M, [1, 2], '<<', opts{:});
    b = -log(I.svd(1));
    x = zeros(length(I.svd));
    y = zeros(length(I.svd));
    for i = 1:length(x)
        x(i) = 2 * sqrt(b * log(I.svd(1) / I.svd(i)));
        y(i) = i - 1;
    end
    toc;
    plot(x, y);
    xlabel('2$\sqrt{bln(\frac{\lambda_{max}}{\lambda})}$', 'Interpreter', 'latex');
    ylabel('n($\lambda$)', 'Interpreter', 'latex');
    savefig('RDMSpec.fig');
end
    