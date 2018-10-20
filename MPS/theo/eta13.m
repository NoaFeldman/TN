function res = eta13(l, t, epsilon)
    % (90) in https://arxiv.org/pdf/1501.00568.pdf
    res = l^2 * epsilon^2 ./ (t.^2 - l^2).^2;
    res(t <= l) = 1;
end