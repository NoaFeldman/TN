function [l, r, I] = myOrthoQS(M, idx, dir, opts)
    if (length(M.Q) == 3)
        [r, l, I] = orthoQS(M, idx, dir, opts{:});
    else
        try
            [l, r, I] = orthoQS(M, idx, dir, opts{:});
        catch E
            disp(E)
            % For what seems to be a very specific matrix in
            % getGroundState(256, -1, 3) there is an error decomposing one
            % block which has quite small numbers anyway. This is some patch
            % that's supposed to overcome it.
            r = QSpace();
            l = QSpace();
            I.svd = [];
            I.svd2tr = 0;
            rightQ = zeros(length(M.Q{1}), 1);
            for i = idx
                rigthQ = rightQ + M.Q{i};
            end
            for q = rightQ
                m = projectSZ(M, q);
                [lm, rm, Im] = orthoQS(m, idx, dir, opts{:});
                r = r + rm;
                l = l + lm;
                I.svd = [I.svd; Im.svd];
                I.svd2tr = I.svd2tr + Im.svd2tr;
            end
            I.svd = sort(I.svd, 'descend');
        end
    end
end