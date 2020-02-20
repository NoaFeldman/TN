F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
nUp = contract(F(1)', '23', F(1), '13');
nDown = contract(F(2)', '23', F(2), '13');
nFull = contract(nUp, 2, nDown, 1);
nUpOnly = cutQSpaceRows(nUp, 1);
nDownOnly = cutQSpaceRows(nDown, 1);
load('/home/noa/TN/MPS/kondoRF/params.mat')