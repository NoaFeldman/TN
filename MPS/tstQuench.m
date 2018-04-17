function tstQuench()
    quenchRDMSpectrum(32, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q32200', 1e-200);
    quenchRDMSpectrum(32, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q3225', 1e-25);
    quenchRDMSpectrum(64, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q64200', 1e-200);
    quenchRDMSpectrum(64, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q6425', 1e-25);
    quenchRDMSpectrum(128, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q128200', 1e-200);
    quenchRDMSpectrum(128, 0, 1, 0, 0, 1e-2, 10, 0, 200, 'Q12825', 1e-25);
end
