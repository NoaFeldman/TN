function [HLR, HLR2] = getHLR(H, psi, l, dir, HLR, HLR2)
    % Returns <H> for the lth site:
    % If l is in the begining of the chain, returns
    %  _
    % | |--
    % | |
    % | |-- for dir='>>', and the miror QSpace for dir = '<<'
    % | |
    % |_|--
    %  
    % else, performs
    %  _        _
    % | |--  --| |--
    % | |      | |  
    % | |--  --| |--
    % | |      | |
    % |_|--  --|_|--
    %
    % HLR.identityChain is I x I x I ... for all sites contracted (two degree tensor)
    % HLR.opSum  is H(1).single x I x I ... + I x H(2).single x I ... +
    %   H.l2r(1) x H.r2l(2) x I ...  (two degree tensor)
    % HLR.openOp is I x I x ... x H(l).l2r  (three degree tensor)
    if (strcmp(dir, '>>'))
        if (l == 0)
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
            HLR2.toClose = QSpace;
            HLR2.toContinue = QSpace;
            return;
        end
        % propagate HLR.opSum by one site, with Identity acting on the new
        % site.
        %  __       _
        % |  |--1--(_)--
        % |  |      |
        % |  |      3
        % |  |      |
        % |__|--2--(_)--
        HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 1), '12', psi(l), '12*');
        % Add the donation of the double operator to opSum, by contracting
        % HLR.openOp with H.r2l(l).
        HLR.opSum  = HLR.opSum + contract( ...
            contract(contract(HLR.openOp, 2, psi(l), 1), 2, psi(l), '1*'), '124', ...
            H.r2l(l), '321');
        % psiPsi = 
        %     ___ 
        %  --(___)--
        % |    |
        % |
        % |   _|_
        %  --(___)--
        psiPsi = contract(psi(l), 1, psi(l), '1*');
        HLR.opSum  = HLR.opSum + contract(H.single(l), '21', psiPsi, '13');
        HLR.opSum = HLR.opSum + contract( ...
            contract(contract(HLR2.toClose, 2, psi(l), 1), 2, psi(l), '1*'), '124', ...
            H.r2l2(l), '321');
        HLR.openOp = contract(H.l2r(l), '21', psiPsi, '13');
        if l > 1
            HLR2.toClose  = contract(contract(HLR2.toContinue, 2, psi(l), 1), '23', psi(l), '12*');
            HLR2.toContinue = contract(H.l2r2(l), '21', psiPsi, '13');
        end
    end
    if (strcmp(dir, '<<'))
        if (l == length(psi)+1)
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
            HLR2.toClose = QSpace;
            HLR2.toContinue = QSpace;
            return;
        end
        % propagate HLR.opSum by one site, with Identity acting on the new
        % site.
        %     _       __
        %  --(_)--1--|  |
        %     |      |  |
        %     3      |  |
        %     |      |  |
        %  --(_)--2--|__|
        HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 3), '13', psi(l), '32*');
        % Add the donation of the double operator to opSum, by contracting
        % HLR.openOp with H.l2r(l).
        HLR.opSum  = HLR.opSum + contract( ...
            contract(contract(HLR.openOp, 2, psi(l), 3), 2, psi(l), '3*'), '135', ...
            H.l2r(l), '321');
        % psiPsi = 
        %     ___ 
        %  --(___)--
        %      |    |
        %           |
        %     _|_   |
        %  --(___)--
        psiPsi = contract(psi(l), 3, psi(l), '3*');
        HLR.opSum  = HLR.opSum + contract(H.single(l), '21', psiPsi, '24');
        HLR.openOp = contract(H.r2l(l), '21', psiPsi, '24');
        
        HLR.opSum = HLR.opSum + contract( ...
            contract(contract(HLR2.toClose, 1, psi(l), 3), 1, psi(l), '3*'), '135', ...
            H.l2r2(l), '321');

        if (l < length(psi)-1)
            HLR2.toClose  = contract(psi(l), '23*', contract(HLR2.toContinue, 1, psi(l), 3), '41', [3 1 2]);
            HLR2.toContinue = contract(psiPsi, '24', H.r2l2(l), '21');        
        end
    end
            
    
    