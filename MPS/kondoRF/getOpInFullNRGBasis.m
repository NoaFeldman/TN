function O = getOpInFullNRGBasis(opKKArr, k, NRGK, HFull)
    OKT = contract(contract(opKKArr(k-1), 2, NRGK.AK, 1), '13', ...
                        NRGK.AT, '13*');
    OTK = contract(contract(opKKArr(k-1), 2, NRGK.AT, 1), '13', ...
                        NRGK.AK, '13*');
    OTT = contract(contract(opKKArr(k-1), 2, NRGK.AT, 1), '13', ...
                        NRGK.AT, '13*');
    OKK = opKKArr(k);
    O = getFullOp(OKK, OKT, OTK, OTT, HFull);
end

function O = getFullOp(OKK, OKT, OTK, OTT, HFull)
    O = 1e-99 * HFull;
    for q = 1:length(O.Q{1})
        kk = find((OKK.Q{1}(:, 1) == O.Q{1}(q, 1)) & ...
                  (OKK.Q{1}(:, 2) == O.Q{1}(q, 2)));
        kt = find((OKT.Q{1}(:, 1) == O.Q{1}(q, 1)) & ...
                  (OKT.Q{1}(:, 2) == O.Q{1}(q, 2)));
        tk = find((OTK.Q{1}(:, 1) == O.Q{1}(q, 1)) & ...
                  (OTK.Q{1}(:, 2) == O.Q{1}(q, 2)));
        tt = find((OTT.Q{1}(:, 1) == O.Q{1}(q, 1)) & ...
                  (OTT.Q{1}(:, 2) == O.Q{1}(q, 2)));
        kLength = 0;
        if ~isempty(kk)
            kLength = length(OKK.data{kk});
            O.data{q}(1:kLength, 1:kLength) = OKK.data{kk};
        end
        if ~isempty(tt)
            tLength = length(OTT.data{tt});
            O.data{q}(kLength + 1:kLength + tLength, kLength + 1:kLength + tLength) = ...
                OTT.data{tt};
        end
        if ~isempty(kk) && ~isempty(tt)
            O.data{q}(kLength + 1:kLength + tLength, 1:kLength) = ...
                OTK.data{tk};
            O.data{q}(1:kLength, kLength + 1:kLength + tLength) = ...
                OKT.data{kt};
        end
    end
end