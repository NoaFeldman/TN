function saveNegs(threadNum, ls, mode, filename)
    maxNumCompThreads(threadNum);
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    load('gs512');
    saveNegativitySpectrum(gs512, ls, mode, filename);
end
