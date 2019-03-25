% How to read [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin',
% 'NC', 1) output: 
% [  0      1     ; 0     -1   ;  0   2   ];
%    N_out  Sz_out; N_in  Sz_in;  N_m Sz_m];

% N takes values of  -1(no particles) , 0(one particle), 1(two particles)
% Sz takes values of -1(one down particle), 0(0/2 particles), 1(one up
% particle)
% This means that:
% F(1) is a_up
% F(2) is a_down
% S(1) is S_+
% S(2) is S_-
% S(3) is S_z
% NUp = contract(F(1)', '23', F(1), '13');
% NDown = contract(F(2)', '13', F(2), '23', [2 1]);

% vac = IS.E;
% vac.data = vac.data(1);
% vac.Q{1} = vac.Q{1}(1, :);
% vac.Q{2} = vac.Q{2}(1, :);
