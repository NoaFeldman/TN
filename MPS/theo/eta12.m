function res = eta12(l, t, epsilon)
    % (89),(91) in https://arxiv.org/pdf/1501.00568.pdf
    res = zeros(1, length(t));
    inds = find(t <=  l);
    res(inds) = (l - t(inds)) ./ (l + t(inds));
    inds = find(t > l);
    res(inds) = l^2 * epsilon^2 ./ (4 .* t(inds).^2 .* (t(inds).^2 - l^2));
end