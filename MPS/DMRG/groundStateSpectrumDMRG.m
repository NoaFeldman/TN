function groundStateSpectrumDMRG(L, JZ, LA)
    if (nargin == 2)
        LA = L/2;
    end
    [psi, ~, ~, ~] = getGroundState(L, 0, 1, JZ, 0);
    saveRDMSpectrum(strcat('groundStateSpectrumDMRG_', int2str(L), '_', int2str(LA), '_', num2str(JZ)), psi, LA);
end