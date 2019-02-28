function overlap = getOverlap(psia, psib)
    % Calculates < psia | psib >
    % k is the index of the site in mixed canonical form.
    if length(psia(1).info.itags) == 3
        curr = contract(psia(1), '12', psib(1), '12*');
    elseif length(psia(1).info.itags) == 4
        curr = contract(psia(1), '123', psib(1), '123*');
    end
    for i = 2:length(psia) - 1
        curr = contract(curr, 1, psia(i), 1);
        if length(psia(1).info.itags) == 3
            curr = contract(curr, '12', psib(i), '12*');
        elseif length(psia(1).info.itags) == 4
            curr = contract(curr, '123', psib(i), '123*');
        end
    end
    curr = contract(curr, 1, psia(length(psia)), 1);
    if length(psia(1).info.itags) == 3
        curr = contract(curr, '123', psib(length(psia)), '123*');
    elseif length(psia(1).info.itags) == 4
        curr = contract(curr, '1234', psib(length(psia)), '1234*');
    end
    overlap = getscalar(curr);
end