function A=ctranspose(A)
% overloading the ' operator
% accept hyper index
% Wb,Nov30,09

  A = builtin('ctranspose', A);

  for k=1:numel(A), Ak=A(k); r=numel(Ak.Q); cgflag=gotCGS(Ak);

     if ~r && isempty(Ak.data) && ~cgflag, continue; end

     if mod(r,2)==0, p=[r/2+1:r, 1:r/2];
     elseif r==3 && isequal(Ak.info.otype,'operator'), p=[2 1 3];
     else error('Wb:ERR','got rank-%d QSpace',r); end
     pc=[sprintf('%g',p) '*'];

   % matlab/2016a requires QSpace() here
     A(k)=QSpace(permuteQS(Ak,pc));
  end

end

