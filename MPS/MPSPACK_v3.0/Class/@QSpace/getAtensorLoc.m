function A=getAtensorLoc(varargin)
% function A=getAtensorLoc(X [,qvac][,l])
%    
%    get Atensor by extracting local state space from local
%    operator X together with vacuum state on virtual bond
%    (default for qvac; optionally any other bond-space can
%    be specified as third argument, instead).
%
%    l \in {'L','R'} specifies the location of vacuum state
%    using LRs index order (default: l='L').
%
%    All arguments except for the last (l) are directly
%    handed over to getIdentityQS(); hence this also allows
%    more general QSpace tensors X and qvac together with
%    the explicit specification of the space index to use.
%
% Wb,Dec28,14

% formerly addsingleton() // Wb,Dec28,14

  if nargin<1 || nargin>5 || numel(varargin{1})~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

% X=varargin{1};
  l='L';

  if nargin>1, q=varargin{end};
     if ischar(q) && (isequal(q,'L') || isequal(q,'R'))
        l=q; varargin=varargin(1:end-1);
     end
  end

  if numel(varargin)>1 && isempty(varargin{end}) % allow to specify qvac=[]
     varargin=varargin(1:end-1);
  end

  n=numel(varargin);
  if n<2 || (n==2 && ~mpsIsQSpace(varargin{2}))
     varargin=[varargin, {getvac(varargin{1},0)}]; % {X [,ix], q }
  end

  A=getIdentityQS(varargin{:}); % index order: loc,vac,combined
     if l=='L', p=[2 3 1]; % -> vac,combined,loc
     else       p=[3 2 1]; % -> combined,vac,loc
     end
  A=permuteQS(A,p);
  A.info.otype='A-matrix';

% A=QSpace(A);
  A=class(A,'QSpace');

end

