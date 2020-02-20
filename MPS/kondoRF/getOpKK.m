function op = getOpKK(NRG, op0)
    % Get an array of impurity operators in each shell. 
    % op(k) = 
    %   __    __    __
    %  |A0|--|op|--|A0|
    %    |           |
    %  |A1|--------|A1|
    %    |           |
    %   ...         ...
    %    |           |
    %  |AK|--------|AK|
    %    |           |
    %
    k = 1;
    op = op0;
    while k <= length(NRG) - 1
        if k == 1
            op(k) = contract(NRG(k).AK, '13*', contract(NRG(k).AK, 3, op0, 2), '13');
        else
            op(k) = contract(NRG(k).AK, '13*', contract(op(k-1), 2, NRG(k).AK, 1), '13');
            op(k) = op(k) + 1e-99 * NRG(k).HK;
        end
        k = k + 1;
    end
    op(k) = contract(NRG(k).AT, '13*', contract(op(k-1), 2, NRG(k).AT, 1), '13');
    op(k) = op(k) + 1e-99 * NRG(k).HT;
end