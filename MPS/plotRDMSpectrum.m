function plotRDMSpectrum(N, Delta1, Delta2, Delta3, Delta4)
    tic;
    Delta = [Delta1, Delta2, Delta3, Delta4];
    mMax = 8;
    weights = zeros(length(Delta), mMax + 1);
    statesNums = zeros(length(Delta), mMax + 1);
    % Recreate Fig.2 from https://arxiv.org/pdf/1407.3779.pdf
    % Added data for Fig.4. If this doesn't work, try: 
    % Verifying l of A, 
    for d = 1:length(Delta)
        [psi, H, HR, HL] = myStartup(N, 0, 1, Delta(d), 0);
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
        eigenstatesCounter = 0;
        while(i > floor(length(psi)/2) + 1)    
            [HL, HR, psi, E, i, M] = dmrgStep(HL, HR, H, psi, i, '<<', opts);
        end
        for m = 0:mMax
            % Sweep to mid chain 
            i = length(psi);
            % force S^z(A) = m
            projector = getProjector(psi(length(psi)/2), m);
            temp = psi(length(psi)/2);
            psi(length(psi)/2) = contract(psi(length(psi)/2), 3, projector, 1);
            if (isempty(psi(length(psi)/2)))
                psi(length(psi)/2) = temp;
                break;
            end
            [~, ~, ~, ~, i, M] = dmrgStep(HL, HR, H, psi, i, '<<', opts);
            [~, ~, I] = orthoQS(M, [1, 2], '<<', opts{:});
            b = -log(I.svd(1));
            weight = 0;
            for i = 1:length(I.svd)
                if (I.svd(i) < 1e-9)
                    break;
                end
                eigenstatesCounter = eigenstatesCounter + 1;
                weight = weight + I.svd(i);
                x(eigenstatesCounter) = 2 * sqrt(b * log(I.svd(1) / I.svd(i)));
            end
            statesNums(d, m+1) = eigenstatesCounter;
            weights(d, m+1) = weight;
            psi(length(psi)/2) = temp;
        end
%        disp(strcat('Finished calculating GS and spectrum for lambda = ', num2str(Delta(d))));
        toc;
    end
%     xlabel('2$\sqrt{bln(\frac{\lambda_{max}}{\lambda})}$', 'Interpreter', 'latex');
%     ylabel('n($\lambda$)', 'Interpreter', 'latex');
%     for i= 1 : 1e3
%         clx(i) = i / 100;
%         cly(i) = besseli(0, (clx(i) / 2)^2);
%         legendInfo{length(Delta) + 1} = ['CL'];
%     end
%     plot(clx, cly);
%     legend(legendInfo);
%     set(gca, 'YScale', 'log');
%     savefig('RDMSpec.fig');
%     hold off;
    save('fig4data', 'weights', 'statesNums', 'x');
end
    
function projector = getProjector(mps, m)
    %
    id = getIdentity(mps, 3);
    mInd = 0;
    for i = 1 : length(id.Q{1})
        if (id.Q{1}(i) == m*2)
            mInd = i;
            break;
        end
    end
    if (mInd == 0)
        projector = QSpace;
        return;
    end
    projector = id;
    for i = 1 : length(projector.Q)
        projector.Q{i} = projector.Q{i}(mInd);
    end
    projector.data = projector.data(mInd);
end