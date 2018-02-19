
% new symmetry names in getLocalSpace
% Pcharge (charge parity, i.e. instead of Acharge)
% Z2charge (same as Pcharge, but with qlabels {0,1}, instead of {1,-1}
% Z2charge can be generalized to ZNcharge with N=2,3,4,...
% Wb,Aug30,16

  [FF,Z,IS]=getLocalSpace('FermionS','Pcharge,SU2spin','NC',1,'-v');
  [FF,Z,IS]=getLocalSpace('FermionS','Z2charge,SU2spin','NC',1,'-v');

  A2=getIdentity(IS.E,IS.E);

