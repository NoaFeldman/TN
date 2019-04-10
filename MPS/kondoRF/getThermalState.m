function [Rho0, EGS] = getThermalState(H, T)
    % This function was written by Matan for H.data already diag()ed
    for i = 1:length(H.data)
        H.data{i} = diag(H.data{i});
    end
    
    Rho0 = H;
    [EGS, bGS] = min(cellfun(@(block) min(diag(block)), H.data)); % find ground state energy and sector
    if T == 0 % ground state
        Rho0.data = {1*diag(diag(H.data{bGS}) == EGS)};
        Rho0.Q = cellfun(@(Q) Q(bGS,:), H.Q,'UniformOutput',0);
    else % thermal state
        Rho0.data = cellfun(@(block) diag(exp(-(diag(block)-EGS)/T)), H.data,'UniformOutput',0);
        Rho0 = Rho0 / trace(Rho0); % normalize
    end
end