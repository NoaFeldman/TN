function A=xitags(A)
% function A=xitags(A)
%
%    skip itags of given QSpace(s)
%
% Wb,Apr14,16

% tags: skipitags, skip_itags

  if nargin>1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  for k=1:numel(A)
     A(k).info.itags=regexprep(A(k).info.itags,'[^\*]*','');
  end

end

