function M = setdiag(M,v,k)
% function M = setdiag(M,v[,k])
%
%     Set diagonal k of matrix M to value (or vector) v (default: k=0)
%     If no output is specified, M is set directly caller.
%      
% Wb,Jul25,05 ; Wb,Aug16,16

  if nargin<2 || nargin>3 
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  if nargin<3, k=0; end

  if ~isnumeric(M) || ndims(M)>2
     error('Wb:ERR','\n   ERR invalid usage'); end
  s=size(M);

% k>0 takes upper diagonals (i.e. same convention as in diag(M,k))
  i=[ max(1,1-k), min(s(1),s(2)-k)]; i=i(1):i(2); % row index (1-based)
  n=numel(i);
  if k<=0, j=0:n-1; else j=k:k+n-1; end           % col index (0-based)

% linearized index
  ii=i+s(1)*j;

  nv=numel(v);
  if nv~=1 && nv~=n || ~isvector(v), error('Wb:ERR',...
    '\n   ERR invalid v (size mismatch %g/%g) !?',nv,n); end

  M(ii)=v;

  if ~nargout, nm=inputname(1);
     if ~isempty(nm)
      % analogous to swap() routine from C.Denham
        assignin('caller',nm,M);
        clear M
     end
  end

end

