function dd=datasize(A)
% Function [D,DD]=datasize(A)
% get full listing of size(data{})
% Wb,Apr24,08

  if nargin~=1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  n=length(A.data); r=length(A.Q); dd=ones(n,r);
  for i=1:n
     s=size(A.data{i});
     dd(i,1:length(s))=s;
  end

end

