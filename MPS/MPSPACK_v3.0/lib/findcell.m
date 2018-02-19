function ii=findcell(c,a)
% Function: ii=findcell(c,a)
%
%    Searches cell c for entries a.
%
% Wb,Oct23,06

  s=size(c); n=prod(s); ii=[];
  for i=1:n
     if isequal(c{i},a), ii=[ii,i]; end
  end

end

