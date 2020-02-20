function blocks = findMatrixBlocks(M)
    skip = 0;
    blocks = {};
    N = length(M);
    for i = 1:N
        for b = 1:length(blocks)
            if ~isempty(find(blocks{b} == i)), skip = 1; end
        end
        if skip == 0
            blocks{end+1} = find(M(i, :) ~=0);
            blocks{end} = union(blocks{end}, find(M(:, i) ~=0));
            j = 1;
            while j <= length(blocks{end})
                blocks{end} = union(blocks{end}, find(M(:, blocks{end}(j)) ~=0));
                blocks{end} = union(blocks{end}, find(M(blocks{end}(j), :) ~=0));
                j = j+1;
            end
        end
        skip = 0;
    end
end