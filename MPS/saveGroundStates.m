function saveGroundStates(Ls)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    for l = Ls
        res.(strcat('gs', int2str(l))) = getGroundState(l, 0, 1, 0, 0);
    end
    save('groundStatesDelta0');
end
        