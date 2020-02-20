function [Rho, EGS] = getThermalState(H, T)    
    Rho = H;
    [EGS, bGS] = min(cellfun(@(block) min(diag(block)), H.data)); % find ground state energy and sector
    if T == 0 % ground state
        Rho.data = {1*diag(diag(H.data{bGS}) == EGS)};
        Rho.Q = cellfun(@(Q) Q(bGS,:), H.Q,'UniformOutput',0);
    else % thermal state
        Rho.data = cellfun(@(block) diag(exp(-(diag(block)-EGS)/T)), H.data,'UniformOutput',0);
        Rho = Rho / trace(Rho); % normalize
    end
end