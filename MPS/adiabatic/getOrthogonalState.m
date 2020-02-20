function res = getOrthogonalState(psi)
    % Maybe this will always work?
    res = psi;
    res(1).data{1} = -psi(1).data{1};
end
    
    