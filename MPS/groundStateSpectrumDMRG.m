function groundStateSpectrumDMRG(L, LA)
    if (nargin == 1)
        LA = L/2;
    end
    [psi, ~, ~, ~] = getGroundState(L, 0, 1, 0, 0);
    saveRDMSpectrum(strcat('groundStateSpectrumDMRG_', int2str(L), '_', int2str(LA)), psi, LA);
end