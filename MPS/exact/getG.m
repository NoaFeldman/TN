function G = getG(sign, cicj, u1, v1, v2)
    % Get G_\pm from (B6) in https://arxiv.org/pdf/1804.00632.pdf
    % Assuming A1 and A2 are adjacent.
    u2 = v1 + 1;
    a = ones(v1 - u1 +1, v1 - u1 +1);
    b = ones(v1 - u1 + 1, v2 - u2 + 1) * sign * 1i;
    c = ones(v2 - u2 + 1, v1 - u1 + 1) * sign * 1i;
    d = ones(v2 - u2 + 1, v2 - u2 + 1);
    G = [a b; c d] .* cicj(u1:v2, u1:v2);
end