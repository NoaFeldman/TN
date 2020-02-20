function [UT, utMatrix] = getTimeReversealUnitary()
    % Get the unitary part of the time reversal operator, e^{-i\pi J_y}.
    y1Eigenvector = [-1 -1i * sqrt(2) 1]/2;
    y1Projector = y1Eigenvector' * y1Eigenvector;
    y0Eigenvector = [1 0 1] / sqrt(2);
    y0Projector = y0Eigenvector' * y0Eigenvector;
    ym1Eigenvector = [-1 1i*sqrt(2) 1]/2;
    ym1Projector = ym1Eigenvector' * ym1Eigenvector;
    utMatrix = exp(-1i * pi * 1) * y1Projector + ...
        exp(-1i * pi * 0) * y0Projector + ...
        exp(-1i * pi * (-1)) * ym1Projector;
    for i = 1:3
        for j = 1:3
            if abs(utMatrix(i, j)) < 1e-15
                utMatrix(i, j) = 0;
            end 
        end
    end
    UT = QSpace();
    UT.info.qtype = '';
    UT.info.otype = '';
%     UT.info.itags = {'', '*', '*'};
    UT.info.itags = {'', '*'};
%     UT.Q = {[], [], []};
    UT.Q = {[], []};
    UT.data = {};
    for inSpin = -1:1
        for outSpin = -1:1
            index = length(UT.Q{1}) + 1;
            UT.Q{1}(index, 1) = 2 * inSpin;
            UT.Q{2}(index, 1) = 2 * outSpin;
%             UT.Q{3}(index, 1) = 2 * (inSpin - outSpin);
            UT.data{index, 1} = [utMatrix(inSpin + 2, outSpin + 2)];
        end
    end
end