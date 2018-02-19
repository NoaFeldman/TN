function A=getId(A,qs,k)
% function A=getId(A,qs [,k])
%
%    reduce QSpace A to rank-2 Id tensor which only has
%    vacuum state / scalar representation with each index.
%    if k is specified, itag{k} will be used.
%
% Wb,Oct18,14

% adapated from QSpace/getvac.m // Wb,Oct02,15

  if isempty(A), error('Wb:ERR','\n   ERR got empty QSpace'); end
  if numel(A)~=1, error('Wb:ERR','\n   ERR single QSpace required'); end

  A=getsub(A,1); r0=numel(A.Q);
  A.data{1}=1;

  qs=reshape(qs,1,[]); l0=0;

  A.Q=A.Q(1:2);
  for i=1:2, A.Q{i}(:)=0; end

  for i=1:size(A.info.cgr,2), q=A.info.cgr(i);
     if ~isempty(q.type) && ~isempty(q.qset) && ~isempty(q.qdir)
     % INIT_SCALAR (see also QSpace/appendScalarSymmetry)
       c=emptystruct(q); % Wb,Dec08,14
       c.type=q.type;

       l=length(q.qset)/r0; l0=l0+l;
       c.qset=zeros(1,2*l);
       if l==numel(qs)
          c.qset=[qs, qs]; j=l0-l+1:l0;
          A.Q{1}(1,j)=qs;
          A.Q{2}(1,j)=qs;
       end

       c.qdir='+-';
       c.cgw=['1' 0];

     A.info.cgr(i)=c; end
  end

  A.info.otype='';

  t=regexprep(A.info.itags,'*','');
  if nargin<3, k=0; % n=0; 
   % picking any existing itag is totally arbitary! // Wb,Jul27,16
   % for i=1:length(t);
   %    l=length(t{i}); if n<l, n=l; k=i; end
   % end
  elseif k>numel(t)
   % NB! k=0 as input implies empty/default itags below
     error('Wb:ERR','\n   ERR index out of bounds (k=%d/%d)',k,numel(t));
  end

  if k, A.info.itags={t{k}, [t{k} '*']};
  else  A.info.itags={'', '*'}; end

end

