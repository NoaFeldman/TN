function Mnew = addrc (M)
% addrc - add one row / column (NAN value)
%         such that tools like surface() plot the entire matrix
% Usage:  Mnew = addrc (M)
%
% Wb,Apr24,02

  s = size(M);

% Mnew = [ M;    zeros(1,  s(2)) ];
% Mnew = [ Mnew, zeros(s(1)+1,1) ];

% Wb,Nov30,03 - values = 0 zero may distort the picture! -> copy last row / col
  Mnew = [ M;    M(end,:)];
  Mnew = [ Mnew, Mnew(:,end) ];

return

