function i=isscalar(A,lflag)
% function i=isscalar(A [,opts])
%
%    Checks whether QSpace is a scalar, meaning rank-0
%    with a single number in data [e.g., a scalar is obtained
%    by tracing/contracting all indices resulting in a single
%    number].
%
% Options
%
%    '-l'  lenient - also consider a QSpace with non-empty Q
%          a scalar, as long as data contains a single number.
%
% Wb,Aug04,08

  nq=numel(A.Q); if nq, nq=size(A.Q{1},1); end
  nd=numel(A.data);

     if nq && nq~=nd, error('Wb:ERR',... % safeguard
        '\n   ERR severe QSpace inconsistency (%g/%g) !?',nq,nd);
     end

  if nd==1, nd=numel(A.data{1}); end
  if nd==1, i=1; else i=0; end
  if ~i || (i && ~nq), return; end

% got i==1 && nq==1
  if nargin>1
     if ~isequal(lflag,'-l'), error('Wb:ERR','invalid usage'); end
  else i=0; end

end

