function Q=contract (A, B, ma, mb, permQ);
% Function - contract
% Usage: Q = contract (A, B, ma, mb, permQ);
%
%   general routine to contract certain
%   set of indizes of two tensors A and B
%
%   Usage 2: generalized trace: contract(A, [i1 i2; i3 i4; ...])
%   traces out i1 with i2, i3 with i4, ...
%
% AWb, F.Verstraete

% check generalized trace: contract(A, [i1,i2], [i3,i4], ...)
  if nargin==2
     if ~isequal(size(A),size(B))
     eval(['help ' mfilename]); error('Wb:ERR',''); end

   % Q=wbtrace(A,B); return % see contract_071214.m
     Q=reshape(A,1,[])*reshape(B,[],1);
     return
  end

  if nargin<4, eval(['help ' mfilename]); return; end

  if iscell(ma), ra=ma{2}; ma=ma{1}; else ra=max(ma); end
  if iscell(mb), rb=mb{2}; mb=mb{1}; else rb=max(mb); end

% add singleton if necessary (e.g. if A or B are scalars)
  sa=size(A); sa(end+1:ra)=1; ra=length(sa);
  sb=size(B); sb(end+1:rb)=1; rb=length(sb);
  
  Ia=1:ra; Ia(ma)=[]; % dim index w/out ma
  Ib=1:rb; Ib(mb)=[];

  AA = permute(A,[Ia ma]);   % reorder: ma indizes last
  BB = permute(B,[Ib mb]);

  if ~isequal(sa(ma),sb(mb))
   % Wb,May16,08: safeguard also for case prod(sa(ma))==prod(sb(mb))
     error('Wb:ERR','misfit in contracted index range [%s; %s]',...
     vec2str(sa(ma)), vec2str(sb(mb)));
  end

% group indizes {Ia}, {ma} and rewrite as regular matrix for multiplication
  Q = reshape(AA,[prod(sa(Ia)) prod(sa(ma))]) * ...
      reshape(BB,[prod(sb(Ib)) prod(sb(mb))]).';

% return to remaining index structure: A first, B last
  snew=[sa(Ia) sb(Ib)]; if isempty(snew), snew=[1 1]; end
  Q = reshape(Q,snew);

% adopt different permutation if desired
  if nargin==5, if ~isempty(permQ), Q=permute(Q,permQ); end; end

end

