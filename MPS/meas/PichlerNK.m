function VR2 = PichlerNK(mergedCopies, transformed, startInd, endInd)
    % TODO debug
    [S,IS]=getLocalSpace('Spin', 0.5,'-A');
    n1 = contract(S(3), '23', S(3)', '13')*2;
    n0 = contract(S(3)', '23', S(3), '13')*2;
    
    % Copy 1
    permuteInds1 = [1 4 2 3];
    permuteInds2 = [1 2 4 3];
    copy1 = 2;
    copy2 = 3;
    % K1_n_K2_m = mod(N_K1, 4) = n, mod(N_K2, 4) = m
    chain_K1_0_K2_0 = getIdentity(transformed(1), 1);
    chain_K1_0_K2_1 = QSpace();
    chain_K1_0_K2_2 = QSpace();
    chain_K1_0_K2_3 = QSpace();
    chain_K1_1_K2_0 = QSpace();
    chain_K1_1_K2_1 = QSpace();
    chain_K1_1_K2_2 = QSpace();
    chain_K1_1_K2_3 = QSpace();
    chain_K1_2_K2_0 = QSpace();
    chain_K1_2_K2_1 = QSpace();
    chain_K1_2_K2_2 = QSpace();
    chain_K1_2_K2_3 = QSpace();
    chain_K1_3_K2_0 = QSpace();
    chain_K1_3_K2_1 = QSpace();
    chain_K1_3_K2_2 = QSpace();
    chain_K1_3_K2_3 = QSpace();

    for i = 1:startInd - 1
        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0, 2, mergedCopies(i), 1);
        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0, '123', mergedCopies(i), '123*', [2, 1]);
    end
    for i = startInd:endInd
        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0, 2, transformed(i), 1);
        chain_K1_0_K2_1 = contract(chain_K1_0_K2_1, 2, transformed(i), 1);
        chain_K1_0_K2_2 = contract(chain_K1_0_K2_2, 2, transformed(i), 1);
        chain_K1_0_K2_3 = contract(chain_K1_0_K2_3, 2, transformed(i), 1);
        chain_K1_1_K2_0 = contract(chain_K1_1_K2_0, 2, transformed(i), 1);
        chain_K1_1_K2_1 = contract(chain_K1_1_K2_1, 2, transformed(i), 1);
        chain_K1_1_K2_2 = contract(chain_K1_1_K2_2, 2, transformed(i), 1);
        chain_K1_1_K2_3 = contract(chain_K1_1_K2_3, 2, transformed(i), 1);
        chain_K1_2_K2_0 = contract(chain_K1_2_K2_0, 2, transformed(i), 1);
        chain_K1_2_K2_1 = contract(chain_K1_2_K2_1, 2, transformed(i), 1);
        chain_K1_2_K2_2 = contract(chain_K1_2_K2_2, 2, transformed(i), 1);
        chain_K1_2_K2_3 = contract(chain_K1_2_K2_3, 2, transformed(i), 1);
        chain_K1_3_K2_0 = contract(chain_K1_3_K2_0, 2, transformed(i), 1);
        chain_K1_3_K2_1 = contract(chain_K1_3_K2_1, 2, transformed(i), 1);
        chain_K1_3_K2_2 = contract(chain_K1_3_K2_2, 2, transformed(i), 1);
        chain_K1_3_K2_3 = contract(chain_K1_3_K2_3, 2, transformed(i), 1);
        
        chain_K1_0_K2_0_ =  ...
            contract(contract(chain_K1_0_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_0_K2_1_ =  ...
            contract(contract(chain_K1_0_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_0_K2_2_ =  ...
            contract(contract(chain_K1_0_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_0_K2_3_ =  ...
            contract(contract(chain_K1_0_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);    
        
        chain_K1_1_K2_0_ =  ...
            contract(contract(chain_K1_1_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_1_K2_1_ =  ...
            contract(contract(chain_K1_1_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_1_K2_2_ =  ...
            contract(contract(chain_K1_1_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_1_K2_3_ =  ...
            contract(contract(chain_K1_1_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_0_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2); 
                
        chain_K1_2_K2_0_ =  ...
            contract(contract(chain_K1_2_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_2_K2_1_ =  ...
            contract(contract(chain_K1_2_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_2_K2_2_ =  ...
            contract(contract(chain_K1_2_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_2_K2_3_ =  ...
            contract(contract(chain_K1_2_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_1_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2); 
         
        chain_K1_3_K2_0_ =  ...
            contract(contract(chain_K1_3_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_3_K2_1_ =  ...
            contract(contract(chain_K1_3_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_0, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_0, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_3_K2_2_ =  ...
            contract(contract(chain_K1_3_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_1, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_1, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2);
        
        chain_K1_3_K2_3_ =  ...
            contract(contract(chain_K1_3_K2_3, copy1, n0, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_3, copy1, n1, 2, permuteInds1), copy2, ...
                n0, 2, permuteInds2) + ...
            contract(contract(chain_K1_3_K2_2, copy1, n0, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2) + ...
            contract(contract(chain_K1_2_K2_2, copy1, n1, 2, permuteInds1), copy2, ...
                n1, 2, permuteInds2); 
            
        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_0_K2_1 = contract(chain_K1_0_K2_1_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_0_K2_2 = contract(chain_K1_0_K2_2_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_0_K2_3 = contract(chain_K1_0_K2_3_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_1_K2_0 = contract(chain_K1_1_K2_0_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_1_K2_1 = contract(chain_K1_1_K2_1_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_1_K2_2 = contract(chain_K1_1_K2_2_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_1_K2_3 = contract(chain_K1_1_K2_3_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_2_K2_0 = contract(chain_K1_2_K2_0_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_2_K2_1 = contract(chain_K1_2_K2_1_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_2_K2_2 = contract(chain_K1_2_K2_2_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_2_K2_3 = contract(chain_K1_2_K2_3_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_3_K2_0 = contract(chain_K1_3_K2_0_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_3_K2_1 = contract(chain_K1_3_K2_1_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_3_K2_2 = contract(chain_K1_3_K2_2_, '123', transformed(i), '123*', [2, 1]);
        chain_K1_3_K2_3 = contract(chain_K1_3_K2_3_, '123', transformed(i), '123*', [2, 1]);
    end
    for i = endInd+1:length(mergedCopies)
        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0, 2, mergedCopies(i), 1);
        chain_K1_0_K2_1 = contract(chain_K1_0_K2_1, 2, mergedCopies(i), 1);
        chain_K1_0_K2_2 = contract(chain_K1_0_K2_2, 2, mergedCopies(i), 1);
        chain_K1_0_K2_3 = contract(chain_K1_0_K2_3, 2, mergedCopies(i), 1);
        chain_K1_1_K2_0 = contract(chain_K1_1_K2_0, 2, mergedCopies(i), 1);
        chain_K1_1_K2_1 = contract(chain_K1_1_K2_1, 2, mergedCopies(i), 1);
        chain_K1_1_K2_2 = contract(chain_K1_1_K2_2, 2, mergedCopies(i), 1);
        chain_K1_1_K2_3 = contract(chain_K1_1_K2_3, 2, mergedCopies(i), 1);
        chain_K1_2_K2_0 = contract(chain_K1_2_K2_0, 2, mergedCopies(i), 1);
        chain_K1_2_K2_1 = contract(chain_K1_2_K2_1, 2, mergedCopies(i), 1);
        chain_K1_2_K2_2 = contract(chain_K1_2_K2_2, 2, mergedCopies(i), 1);
        chain_K1_2_K2_3 = contract(chain_K1_2_K2_3, 2, mergedCopies(i), 1);
        chain_K1_3_K2_0 = contract(chain_K1_3_K2_0, 2, mergedCopies(i), 1);
        chain_K1_3_K2_1 = contract(chain_K1_3_K2_1, 2, mergedCopies(i), 1);
        chain_K1_3_K2_2 = contract(chain_K1_3_K2_2, 2, mergedCopies(i), 1);
        chain_K1_3_K2_3 = contract(chain_K1_3_K2_3, 2, mergedCopies(i), 1);

        chain_K1_0_K2_0 = contract(chain_K1_0_K2_0, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_0_K2_1 = contract(chain_K1_0_K2_1, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_0_K2_2 = contract(chain_K1_0_K2_2, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_0_K2_3 = contract(chain_K1_0_K2_3, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_1_K2_0 = contract(chain_K1_1_K2_0, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_1_K2_1 = contract(chain_K1_1_K2_1, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_1_K2_2 = contract(chain_K1_1_K2_2, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_1_K2_3 = contract(chain_K1_1_K2_3, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_2_K2_0 = contract(chain_K1_2_K2_0, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_2_K2_1 = contract(chain_K1_2_K2_1, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_2_K2_2 = contract(chain_K1_2_K2_2, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_2_K2_3 = contract(chain_K1_2_K2_3, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_3_K2_0 = contract(chain_K1_3_K2_0, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_3_K2_1 = contract(chain_K1_3_K2_1, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_3_K2_2 = contract(chain_K1_3_K2_2, '123', mergedCopies(i), '123*', [2, 1]);
        chain_K1_3_K2_3 = contract(chain_K1_3_K2_3, '123', mergedCopies(i), '123*', [2, 1]);
    end
    % TODO works without the minus signs... check
    VR2 = abs(getscalar(chain_K1_0_K2_0)) + ...
         -abs(getscalar(chain_K1_0_K2_2)) + ...
          abs(getscalar(chain_K1_1_K2_1)) + ...
         -abs(getscalar(chain_K1_1_K2_3)) + ...
         -abs(getscalar(chain_K1_2_K2_0)) + ...
          abs(getscalar(chain_K1_2_K2_2)) + ...
         -abs(getscalar(chain_K1_3_K2_1)) + ...
          abs(getscalar(chain_K1_3_K2_3));
end