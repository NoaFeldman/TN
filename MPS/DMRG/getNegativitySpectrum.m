function spectrum = getNegativitySpectrum(E12T2)
    E12T2 = contract(E, '13', getIdentity(E12T2, 2, E12T2, 3), '12');
    E12T2 = contract(E12T2, '12', getIdentity(E12T2, 1, E12T2, 2), '12');
    
end