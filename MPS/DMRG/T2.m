function partialTransposed = T2(E12)
    % Partially transpose E12 from Eq 64 in https://arxiv.org/pdf/1605.00674.pdf
    %      __                __
    %  a--|  |--c        b--|  |--c
    %     |  |     >>>>     |  |
    %  b--|__|--d        a--|__|--d
    %
    partialTransposed = E12;
    partialTransposed.Q{1, 1} = E12.Q{1, 2};
    partialTransposed.Q{1, 2} = E12.Q{1, 1};
    tags = E12.info.itags;
    temp = tags{1};
    tags{1} = tags{2};
    tags{2} = temp;
    partialTransposed.info.itags = tags;
end