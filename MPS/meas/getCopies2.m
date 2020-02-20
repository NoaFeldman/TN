function copies = getCopies2(psi)
    copies = QSpace(length(psi));
    for j = 1:length(psi)
        curr = psi(j);
        curr.Q = [curr.Q {zeros(length(curr.Q{1}), 1)}];
        curr.info.itags{4} = '';
        curr2 = curr;
        curr2.info.itags{4} = '*';
        curr = contract(curr, 4, curr2, 4);
        curr.info.itags = {strcat('A', curr.info.itags{1}), ...
                           strcat('A', curr.info.itags{2}), ...
                           strcat('A', curr.info.itags{3}), ...
                           strcat('B', curr.info.itags{1}), ...
                           strcat('B', curr.info.itags{2}), ...
                           strcat('B', curr.info.itags{3})};
        copies(j) = curr;
    end
end
                       
