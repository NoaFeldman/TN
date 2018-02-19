function A=transpose(A)
% overloading the .' operator (element wise transpose)

  if isempty(A) || numel(A)==1 && isempty(A.Q) && isempty(A.data)
     return
  end

  A=builtin('transpose',A);
  % if numel(A)>1 && builtin('ndims',A)==2
  %    A=builtin('permute',A,[2 1]);
  % end

  for k=1:numel(A), r=numel(A(k).Q);
     if r==2, p=[2 1];
     elseif r==3 && isequal(A(k).info.otype,'operator'), p=[2 1 3];
     else error('Wb:ERR',sprintf('got rank-%d QSpace',r)); end
     A(k)=permuteQS(A(k),p);
  end

end

