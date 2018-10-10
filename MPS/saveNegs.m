function saveNegs()
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    load('gs512');
    saveNegativities(gs512, [8 16 32 48 64], 3:5, 'negs512');
end
