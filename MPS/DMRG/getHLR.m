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
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end
    if (strcmp(dir, '>>'))
        if (l == 0)
            HLR.opSum = QSpace;
            if strcmp(bc, 'open')
                HLR.openOp = QSpace;
                HLR2.toClose = QSpace;
                HLR2.toContinue = QSpace;
            else
                HLR.openOpUp = QSpace;
                HLR2.toCloseUp = QSpace;
                HLR2.toContinueUp = QSpace;                
                HLR.openOpDown = QSpace;
                HLR2.toCloseDown = QSpace;
                HLR2.toContinueDown = QSpace;                
            end
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
        if strcmp(bc, 'open')
            HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 1), '12', psi(l), '12*');
            % Add the donation of the double operator to opSum, by contracting
            % HLR.openOp with H.r2l(l).
            HLR.opSum  = HLR.opSum + contract( ...
                contract(contract(HLR.openOp, 2, psi(l), 1), 2, psi(l), '1*'), '124', ...
                H.r2l(l), '321');
        else
            HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 1), '123', psi(l), '123*');
            % Add the donation of the double operator to opSum, by contracting
            % HLR.openOp with H.r2l(l).
            HLR.opSum  = HLR.opSum + contract( ...
                contract(contract(HLR.openOpUp, 2, psi(l), 1), 2, psi(l), '1*'), '12356', ...
                H.r2lUp(l), '52413');
            HLR.opSum  = HLR.opSum + contract( ...
                contract(contract(HLR.openOpDown, 2, psi(l), 1), 2, psi(l), '1*'), '12356', ...
                H.r2lDown(l), '52413');
        end
        % psiPsi = 
        %     ___ 
        %  --(___)--
        % |    |
        % |
        % |   _|_
        %  --(___)--
        psiPsi = contract(psi(l), 1, psi(l), '1*');
        if strcmp(bc, 'open')
            HLR.opSum  = HLR.opSum + contract(H.single(l), '21', psiPsi, '13');            
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toClose, 2, psi(l), 1), 2, psi(l), '1*'), '124', ...
                H.r2l2(l), '321');
        else
            HLR.opSum  = HLR.opSum + contract(H.single(l), '2143', psiPsi, '1425');
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toCloseUp, 2, psi(l), 1), 2, psi(l), '1*'), '12356', ...
                H.r2l2Up(l), '52413');
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toCloseDown, 2, psi(l), 1), 2, psi(l), '1*'), '12356', ...
                H.r2l2Down(l), '52413');

        end
        if strcmp(bc, 'open')
            HLR.openOp = contract(H.l2r(l), '21', psiPsi, '13');
        else
            HLR.openOpUp = contract(H.l2rUp(l), '2143', psiPsi, '1425');
            HLR.openOpDown = contract(H.l2rDown(l), '2143', psiPsi, '1425');
        end
        if strcmp(bc, 'open')        
            if l > 1
                HLR2.toClose  = contract(contract(HLR2.toContinue, 2, psi(l), 1), '23', psi(l), '12*');
            end
            HLR2.toContinue = contract(H.l2r2(l), '21', psiPsi, '13');
        else
            HLR2.toCloseUp  = contract(contract(HLR2.toContinueUp, 2, psi(l), 1), '234', psi(l), '123*');
            HLR2.toContinueUp = contract(H.l2r2Up(l), '2143', psiPsi, '1425');                
            HLR2.toCloseDown  = contract(contract(HLR2.toContinueDown, 2, psi(l), 1), '234', psi(l), '123*');
            HLR2.toContinueDown = contract(H.l2r2Down(l), '2143', psiPsi, '1425');                
        end
    end
    if (strcmp(dir, '<<'))
        if (l == length(psi)+1)
            HLR.opSum = QSpace;
            if strcmp(bc, 'open')
                HLR.openOp = QSpace;
                HLR2.toClose = QSpace;
                HLR2.toContinue = QSpace;
            else
                HLR.openOpUp = QSpace;
                HLR2.toCloseUp = QSpace;
                HLR2.toContinueUp = QSpace;                
                HLR.openOpDown = QSpace;
                HLR2.toCloseDown = QSpace;
                HLR2.toContinueDown = QSpace;                
            end
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
        if strcmp(bc, 'open')
            HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 3), '13', psi(l), '32*');
            % Add the donation of the double operator to opSum, by contracting
            % HLR.openOp with H.l2r(l).
            HLR.opSum  = HLR.opSum + contract( ...
                contract(contract(HLR.openOp, 2, psi(l), 3), 2, psi(l), '3*'), '135', ...
                H.l2r(l), '321');
        else
            HLR.opSum = contract(contract(HLR.opSum, 1, psi(l), 4), '134', psi(l), '423*');
            % Add the donation of the double operator to opSum, by contracting
            % HLR.openOp with H.l2r(l).
            HLR.opSum  = HLR.opSum + contract( ...
                contract(contract(HLR.openOpUp, 2, psi(l), 4), 2, psi(l), '4*'), '13647', ...
                H.l2rUp(l), '52143');
        end
        % psiPsi = 
        %     ___ 
        %  --(___)--
        %      |    |
        %           |
        %     _|_   |
        %  --(___)--
        if strcmp(bc, 'open')
            psiPsi = contract(psi(l), 3, psi(l), '3*');
            HLR.opSum  = HLR.opSum + contract(H.single(l), '21', psiPsi, '24');
            HLR.openOp = contract(H.r2l(l), '21', psiPsi, '24');
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toClose, 2, psi(l), 3), 2, psi(l), '3*'), '135', ...
                H.l2r2(l), '321');
        else
            psiPsi = contract(psi(l), 4, psi(l), '4*');
            HLR.opSum  = HLR.opSum + contract(H.single(l), '2143', psiPsi, '2536');
            HLR.openOpUp = contract(H.r2lUp(l), '2143', psiPsi, '2536');
            HLR.openOpDown = contract(H.r2lDown(l), '2143', psiPsi, '2536');
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toCloseUp, 1, psi(l), 4), 1, psi(l), '4*'), '13647', ...
                H.l2r2Up(l), '52143');
            HLR.opSum = HLR.opSum + contract( ...
                contract(contract(HLR2.toCloseDown, 1, psi(l), 4), 1, psi(l), '4*'), '13647', ...
                H.l2r2Down(l), '52143');
        end

        if strcmp(bc, 'open')
            if (l < length(psi)-1)
                HLR2.toClose  = contract(psi(l), '23*', contract(HLR2.toContinue, 2, psi(l), 3), '42', [2 3 1]);
            end
            HLR2.toContinue = contract(psiPsi, '24', H.r2l2(l), '21', [3 1 2]);        
        else
            HLR2.toCloseUp  = contract(psi(l), '234*', contract(HLR2.toContinueUp, 2, psi(l), 4), '452', [3 1 2]);                
            HLR2.toCloseDown  = contract(psi(l), '234*', contract(HLR2.toContinueDown, 2, psi(l), 4), '452', [3 1 2]);                
            HLR2.toContinueUp = contract(psiPsi, '2536', H.r2l2Up(l), '2143', [3 1 2]);        
            HLR2.toContinueDown = contract(psiPsi, '2536', H.r2l2Down(l), '2143', [3 1 2]);        
        end
    end
end
            
    
    