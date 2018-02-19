function A=getsub(A,varargin)
% function A=getsub(A,idx)
%    select given subspace where idx represents direct index
%    into records of QSpace.
%
% function A=getsub(A,q,dim)
%    select subspace that has given symmetry/ies in given dimension
%
% Wb,Aug20,08 ; Wb,Jun19,13

  if nargin==2, idx=varargin{1};
   % NB! DO NOT use unique() as the idx-order must be preserved
   % eg. to match block order of different QSpaces etc.!
   % ensure vector structure though // Wb,Jun04,16
     if ~isvector(idx)
        wblog('WRN','got non-vector as index (%s) !?',sizestr(idx));
        idx=idx(:);
     end
     if any(diff(sort(idx))==0)
        error('Wb:ERR','\n   ERR invalid usage (got non-unique idx)'); 
     end
  elseif nargin==3
     q=varargin{1}; k=varargin{2};
     [ia,ib,Im]=matchIndex(cat(2,A.Q{k}),q);
     idx=unique(ia);
  else
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  for p=1:length(A.Q), A.Q{p}=A.Q{p}(idx,:); end

% reduce cgr first (since gotCGS checks data length)
  if isfield(A.info,'cgr') && ~isempty(A.info.cgr)
     A.info.cgr=A.info.cgr(idx,:);
  end

  A.data=A.data(idx);

end

