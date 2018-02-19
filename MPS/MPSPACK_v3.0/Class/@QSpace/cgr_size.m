function s=cgr_size(A,i,j)
% function s=cgr_size(A,i,j)
% Wb,Jan08,12

  if nargin~=3 || numel(A)~=1 || numel(i)~=1 || numel(j)~=1
  error('Wb:ERR','\n   ERR invalid usage'); end

  r=numel(A.Q);
  if ~gotCGS(A), s=ones(1,r); return; end

  c=A.info.cgr;
  if iscell(c), c=c{i,j};
  else c=c(i,j); end

  if isfield(c,'size')
     s=c.size;
  else
     if ~isnumeric(c), error('Wb:ERR','\n   ERR invalid cgr data'); end
     s=size(c);
  end

% adjust for trailing singletons
  n=numel(s);
  if n>r
   % NB! may have outer multiplicity
     if n>r+1, error('Wb:ERR',...
       '\n   ERR invalid size (rank too large: %d/%d)',n,r);
     end
     s=s(1:r);
  elseif n<r, s(end+1:r)=1; end
   
end

