% Check that the rAlphas vector contains right eigenvectors of the
% Liouvillian
Lrs = -1i * (contract(LiouB.H, 2, rAlphas, 1) - contract(rAlphas, 2, LiouB.H, 1, [1 3 2])) ...
    + contract(contract(LiouB.L, 2, rAlphas, 1), 2, LiouB.L', 1, [1 3 2]) ...
    - 1/2 * (contract(LiouB.LdagL, 2, rAlphas, 1) + contract(rAlphas, 2, LiouB.LdagL, 1, [1 3 2]));

for alpha = 1:size(rAlphas.data{1}, 3)
    if abs(lambdaAlphas(alpha)) < 1e-12
        if  max(max(abs(Lrs.data{1}(:, :, alpha) ...
                ./Lrs.data{1}(:, :, alpha)))) > 1e-12
           disp(alpha);
           disp(['lambda = ' num2str(lambdaAlphas(alpha))]);
           disp(Lrs.data{1}(:, :, alpha));
        end
    elseif max(max(abs(rAlphas.data{1}(:, :, alpha) * lambdaAlphas(alpha) - ...
               Lrs.data{1}(:, :, alpha)))) > 1e-12
           disp(alpha);
           disp(rAlphas.data{1}(:, :, alpha) / lambdaAlphas(alpha));
           disp(Lrs.data{1}(:, :, alpha));
    end
end

% Check that the lDagAlphas vector contains left eigenvectors of the
% Liouvillian
Lls = 1i * (contract(LiouB.H, 2, lDagAlphas, 1) - contract(lDagAlphas, 2, LiouB.H, 1, [1 3 2])) ...
    + contract(contract(LiouB.L', 2, lDagAlphas, 1), 2, LiouB.L, 1, [1 3 2]) ...
    - 1/2 * (contract(LiouB.LdagL, 2, lDagAlphas, 1) + contract(lDagAlphas, 2, LiouB.LdagL, 1, [1 3 2]));

for alpha = 1:size(rAlphas.data{1}, 3)
    if abs(lambdaAlphas(alpha)) < 1e-12
        if max(max(abs(Lls.data{1}(:, :, alpha) ...
                        ./lDagAlphas.data{1}(:, :, alpha)))) > 1e-12 % This is a little flaky, try to think of a better criterion for matrix == 0.
            disp(alpha);
            disp(['lambda = ' num2str(lambdaAlphas(alpha))]);
            disp(Lls.data{1}(:, :, alpha));
        end
    elseif max(max(abs(lDagAlphas.data{1}(:, :, alpha) * conj(lambdaAlphas(alpha)) - ...
               Lls.data{1}(:, :, alpha)))) > 1e-12
            disp(alpha);
            disp(lDagAlphas.data{1}(:, :, alpha) / lambdaAlphas(alpha));
            disp(Lls.data{1}(:, :, alpha));
    end
end