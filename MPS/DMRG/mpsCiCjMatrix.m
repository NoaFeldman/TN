function [cici, cicicjcj] = mpsCiCjMatrix(psi, S)
    cici = zeros(length(psi), 1);
    cicicjcj = zeros(length(psi));
    N = contract(S(2)', '23', S(2), '13')*2;
    for i = 1:length(psi)
        psiCurr = psi;
        psiCurr(i) = contract(psiCurr(i), 2, N, 2, [1 3 2]);
        cici(i) = expectationVal(psiCurr);
        cicicjcj(i, i) = cici(i);
    end
    for i = 1:length(psi)
        for j = i+1:length(psi)
            psiCurr = psi;
            psiCurr(i) = contract(psiCurr(i), 2, N, 2, [1 3 2]);
            psiCurr(j) = contract(psiCurr(j), 2, N, 2, [1 3 2]);
            cicicjcj(i, j) = expectationVal(psiCurr);
            cicicjcj(j, i) = cicicjcj(i, j);
        end
    end
end

function val = expectationVal(psi)
    chain = contract(psi(1), '12', psi(1), '12*');
    for i = 2:length(psi) - 1
        chain = contract(contract(chain, 1, psi(i), 1), '12', psi(i), '12*');
    end
    val = getscalar(...
        contract(contract(chain, 1, psi(length(psi)), 1), '123',...
                 psi(length(psi)), '123*'));
end