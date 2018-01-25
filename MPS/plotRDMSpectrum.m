function plotRDMSpectrum(N, Delta1, Delta2, Delta3, Delta4)
    tic;
    Delta = [Delta1, Delta2, Delta3, Delta4];
    % Recreate plot fFig.2 from https://arxiv.org/pdf/1407.3779.pdf
    for d = 1:length(Delta)
        [psi, H, HR, HL] = myStartup(N, 0, 1, Delta(d));
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
        for i = 1:length(I.svd)
            x(i) = 2 * sqrt(b * log(I.svd(1) / I.svd(i)));
            y(i) = i - 1;
            if (x(i) > 10)
                break;
            end
        end
        hold on;
        plot(x, y);
        legendInfo{d} = ['Delta = ' num2str(Delta(d))];
        disp(strcat('Finished calculating GS and spectrum for lambda = ', num2str(Delta(d))));
        toc;
    end
    for d = 1:length(Delta)
    end
    xlabel('2$\sqrt{bln(\frac{\lambda_{max}}{\lambda})}$', 'Interpreter', 'latex');
    ylabel('n($\lambda$)', 'Interpreter', 'latex');
    for i= 1 : 1e3
        clx(i) = i / 100;
        cly(i) = besseli(0, (clx(i) / 2)^2);
        legendInfo{length(Delta) + 1} = ['CL'];
    end
    plot(clx, cly);
    legend(legendInfo);
    set(gca, 'YScale', 'log');
    savefig('RDMSpec.fig');
    hold off;
end
    