function saveTDependentNegativities(inDirName, outDirName, l, ns, steps)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    for step : steps
        load(strcat(inDirName, '/psiAtStep', int2str(step)));
        saveNegativities(psi, l, ns, stract(outDirName, '/negsAtStep', int2str(step)));
    end
end