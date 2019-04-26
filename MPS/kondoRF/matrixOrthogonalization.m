function newLAlphas = matrixOrthogonalization(rAlphas, lAlphas, lambdaAlphas)
    % Orthogonalize left and right eigenmatrices based on equation 6 in
    % https://arxiv.org/pdf/1811.03518.pdf.
    newLAlphas = QSpace(length(lAlphas));
    for i = 1:length(rAlphas)
        for j = 1:length(newLAlphas)
            if abs(trace(rAlphas(i)*lAlphas(j))) > 0.5
                if (i ~= j)
                    k = 1;
                end
                newLAlphas(i) = lAlphas(j) / trace(rAlphas(i)*lAlphas(j));
            end
        end
    end
end