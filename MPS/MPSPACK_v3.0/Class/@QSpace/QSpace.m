function [A,varargout] = QSpace(varargin)
% Function [A,varargout] = QSpace(...)
%
%    class constructor for abelian quantum numbers
%
% Usage:
%
%  1: A = QSpace(Q, data [,'check'])
%     A = QSpace(A [,'check']) - optinally checks QSpace structure
%
%  2: A = QSpace(q1, data1, q2, data2, ...) % see initQSpace
%
%  3: A = initQSpace(Q1, M, DB)
%         square operator M with row/col Q1 block size given by DB
%
%  4: A = initQSpace(Q1,Q2,M,'map')
%         matrix M with row/col QIDX given by Q1/Q2
%         eg. used for projectors
%
%  5: A = QSpace(Q, M, 'operator');
%  6: A = QSpace(A, 'unity');
%
%  7: A0=initQSpace(QL,Qs1,Qs2,...); setup identity tensor
%         for product spaces (QL,Qs1,...)
%  
%         it must hold: size(QL,2) = size(Qs1,2) = size(Qs2,2) = ... = QDIM.
%         Qs and QL span the full Hilbert space, hence certain quantum
%         lables may appear several times.
%
% Wb,May12,06  Wb,Aug10,06  Wb,Aug06,07

  if nargin==0, A=QSpace({},{}); return, end
  if nargin==1 && isa(varargin{1},'QSpace')
     A=varargin{1}; return
  end

  if nargin==1 && isempty(varargin{1}), A=QSpace({},{}); return, end

  varargout=cell(1,nargout-1);

  if iscell(varargin{1}) || isstruct(varargin{1}) 
   % usage #1
     switch nargin

      case 1

         if iscell(varargin{1})
            sa=size(varargin{1}); 
            for k=prod(sa):-1:1, A(k)=varargin{1}{k}; end
            A=reshape(A,sa);
         else
            A=varargin{1};
         end

      otherwise
         if nargin==2
          % need extra {} brackets here to get a single structure
            A=struct('Q',{varargin{1}}, 'data',{varargin{2}(:)}, ...
             'info',{[]});
         else
            A=struct('Q',{varargin{1}}, 'data',{varargin{2}(:)}, ...
             'info',{varargin{3}});
         end

         if ~isempty(A.Q) && size(A.Q{1},1)~=length(A.data)
         error('Wb:ERR','QSpace init - size inconsistency'); end
     end

     if ~isfield(A,'info') && numel(A), A(1).info=[]; end

     if ~isequal(varargin{end},'check')
        A=class(A,'QSpace');
        return
     end

  elseif isnumeric(varargin{1})

   % for initialization (like zeros())
     if nargin<=2 && isnumeric(varargin{end});
        if nargin==2 && isscalar(varargin{1}) && isscalar(varargin{2})
            d=cat(2,varargin{:});
        elseif nargin==1 && isvector(varargin{1})
            d=varargin{1};
        else d=[]; end

        if ~isempty(d)
           if length(d)==1, d=[d 1]; end
           A=repmat(QSpace({},{}),d);
           return
        end
     end

   % usage #7
     if isequal(varargin{end},'identity')
        [A,varargout{:}]=initQSpaceA0(varargin{1:end-1});
   % usage #5
     elseif isequal(varargin{end},'operator')
        [A,varargout{:}]=initQSpaceH0(varargin{1:end-1});
   % usage #4
     elseif isequal(varargin{end},'map')
        [A,varargout{:}]=initQMap(varargin{1:end-1});
   % usage #8 QSpace(q1,rank,data,'nosym') -- Wb,Sep20,10
     elseif isequal(varargin{end},'nosym')
      % QSpace(qq,rank,data,'nosym');
        if nargin~=4, error('Wb:ERR','invalid usage #8'); end
        if ~isempty(varargin{3})
            A=struct(...
             'Q',{repmat(varargin(1),1,varargin{2})}, ...
             'data',{varargin(3)}, 'info',{[]});
        else
            A=struct('Q',{{}},'data',{cell(0,1)},'info',{[]});
        end
   % usage #3
     elseif nargin==3 && diff(size(varargin{2}))==0 % i.e. ismatrix
        A=initQSpaceMD(varargin{:});
     else
   % usage #2
        A=initQSpace(varargin{:});
     end

  elseif isa(varargin{1},'QSpace') || isstruct(varargin{1})

   % usage #6
     if isequal(varargin{end},'unity')
        A=initQSpaceUnity(varargin{1:end-1});
     else
        error('Wb:ERR','\n%s ERR invalid QSpace constructor set',lineno);
     end
  else
   % eval(['help ' mfilename]);
     error('Wb:ERR','\n%s ERR invalid QSpace constructor set',lineno);
  end

  for i=1:prod(size(A))
     [isq,s]=mpsIsQSpace(A(i));
     if ~isq, A(i), s, error('Wb:ERR', ...
     sprintf('%s ERR Invalid QSpace constructor set (%d)',lineno,i)); end
  end

  if ~isfield(A,'info') && numel(A), A(1).info=[]; end
  A=class(A,'QSpace');

end % EOM

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function A = initQSpace(varargin)
% initQSpace - initialize QSpace from partial data
%
% Usage: A = initQSpace(q1, data1, q2, data2, ...)
% Example: rank-2 object with two abelian quantum numbers (QDIM=2)
%
%     A=initQSpace(...
%       [-1  0; 0  1],  1, ...
%       [ 0 -1; 1  0], -1, ...
%     );
%
% Wb,Aug08,06

  if nargin<2 || mod(nargin,2)~=0
     eval(['help ' mfilename]);
     error('Wb:ERR', 'Invalid number of arguments (usage #2)');
  end

  for i=1:2:nargin
     if ~isequal(size(varargin(i)), size(varargin(1)))
        eval(['help ' mfilename]);
        error('Wb:ERR','Size mismatch in Q data (usage #2)');
     end
     if ~isnumeric(varargin{i}) || ~isnumeric(varargin{i+1})
        eval(['help ' mfilename]);
        error('Wb:ERR','Invalid data (must be numeric, usage #2)');
     end
  end

  [rank,qdim]=size(varargin{1});
  n=nargin/2;

  qq = permute(reshape( cat(2,varargin{1:2:end}), ...
              [rank,qdim,n]), ...
       [3 2 1]);
  Q = cell(1,rank); for i=1:rank, Q{i}=qq(:,:,i); end

% need extra {} brackets here to get a single structure
  A=struct('Q',{Q}, 'data',{varargin(2:2:end)'}, 'info',{[]});

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function A = initQSpaceMD(varargin)
% initQSpaceMD - initialize QSpace from matrix
%
% Usage: A = QSpace(Q1, M, DB)
%
%    Matrix M must be square with consecutive block
%    dimensions given by DB
%
% Wb,Sep08,06

  if nargin~=3
     eval(['help ' mfilename]);
     error('Wb:ERR', 'Invalid number of arguments (usage #3)');
  end

  Q=varargin{1};
  M=varargin{2};
  D=varargin{3}; n=length(D);

  if size(Q,1)~=n || any(size(M)~=sum(D))
     eval(['help ' mfilename]);
     error('Wb:ERR', 'Dimension mismatch of arguments (usage #3)');
  end

  M = mat2cell(M, D, D);
  mark=zeros(n,n);

  for i=1:n
  for j=1:n, mark(i,j) = any(M{i,j}(:)); end, end

% remember all diagonal blocks, if at least one diagonal block is set
  if any(diag(mark)), mark=mark+eye(size(mark)); end

  [I,J]=find(mark);

% need extra {} brackets here to get a single structure
  A=struct('Q',{{Q(I,:), Q(J,:)}}, 'data', {M(find(mark))}, 'info',{[]});

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [A,K] = initQSpaceH0(varargin)
% initQSpaceH0 - initialize QSpace from matrix
%
% Usage: [A [,K]] = QSpace(Q, M)
%
%    Matrix M must be square with quantum numbers
%    of every index specified in Q
%
% Wb,Sep08,06

  dflag=[];

  if nargin==3
     switch varargin{end}
       case 'diag', dflag=1; varargin=varargin(1:end-1);
       case 'econ', dflag=0; varargin=varargin(1:end-1);
       otherwise e=1;
     end
  end

  if length(varargin)~=2
     error('Wb:ERR', 'invalid operator initialization');
  elseif size(varargin{1},1)~=size(varargin{2},1)
     eval(['help ' mfilename]);
     error('Wb:ERR', 'Dimension mismatch of arguments (usage #4)');
  end

  [Q,K,D]=uniquerows(varargin{1}); n=length(K);
  mark=zeros(n,n); M=cell(n,n);

  for i=1:n, for j=1:n
     a=varargin{2}(K{i},K{j});
     if any(a(:))
        M{i,j}=full(a); mark(i,j)=1;
     end
  end, end

% make sure, index sets in K are vectors sorted in ascending order
  for g=1:length(K)
     p=sort(reshape(K{g},1,prod(size(K{g})))); if ~isequal(p,K{g})
       wblog('WRN','index sets not sorted in uniquerows() ???');
       K{g}=p;
     end
  end

% remember all diagonal blocks, if at least one diagonal block is set
% e.g. relevant for H0 since local basis of A0 will be compared to it
% which better be complete
  if isempty(dflag) && any(diag(mark)) || ~isempty(dflag) && dflag
     for i=1:n
        if mark(i,i), continue; end
        M{i,i}=full(varargin{2}(K{i},K{i}));
        mark(i,i)=1;
     end
  end

  [I,J]=find(mark);

% need extra {} brackets here to get a single structure
  A=struct('Q',{{Q(I,:), Q(J,:)}}, 'data', {M(find(mark))}, 'info',{[]});

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function varargout = initQMap(varargin)
% initQMap - initialize QSpace from matrix
%
% Usage: [A [,K]] = QSpace(Q1,Q2,...,M)
%
%    M is regular matrix
%    with quantum numbers of every index (i,j) specified in Q1,Q2
%
% generalized to arbitrary rank>=2 // Wb,Mar25,14
% Wb,Aug06,07

  n=nargin; if n<3, e=-100;
  else e=0; M=varargin{end};
     for i=1:n
        if ~isnumeric(varargin{i}), e=i;
        elseif i<n && size(varargin{i},1)~=size(M,i), e=100*i; end
     end
  end
  if e
     eval(['help ' mfilename]); error('Wb:ERR', ...
    'dimension mismatch of arguments (initQMap [e=%g])',e);
  end

  n=n-1; s=zeros(1,n);
  for i=1:n
    [Q{i},K{i},D{i}]=uniquerows(varargin{i});
    s(i)=size(Q{i},1);
  % make sure, index sets in K are vectors sorted in ascending order
    for j=1:s(i)
       p=sort(reshape(K{i}{j},1,[])); if ~isequal(p,K{i}{j})
         wblog('WRN','index sets not sorted in uniquerows() ???');
         K{i}{j}=p;
       end
    end
  end

  mark=zeros(s); data=cell(s); i=cell(1,n); I=cell(1,n);
  if n>2, M=full(M); end

  for l=1:prod(s), [i{:}]=ind2sub(s,l);
     for j=1:n, I{j}=K{j}{i{j}}; end
     a=M(I{:});
     if any(a(:))
        data{i{:}}=full(a); mark(i{:})=1;
     end
  end

  l=find(mark); [i{:}]=ind2sub(s,l);
  for j=1:n, I{j}=Q{j}(i{j},:); end

% need extra {} brackets here to get a single structure
  A=struct('Q',{I}, 'data', {data(l)}, 'info',{[]});

  if nargout<2
       varargout={A};
  else varargout={A,K{:}}; end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% build identity operator from other operator basis

function E = initQSpaceUnity(A)

  if nargin~=1 || ~isa(A,'QSpace') && (~isstruct(A) || ~isfield(A,'Q'))
     error('Wb:ERR','\ninvalid usage of %s\n%s', ...
     lineno('%S'), lineno('all'));
  end
  if numel(A)~=1, e=numel(A);
  elseif numel(A.Q)~=2, e=numel(A.Q); else e=0; end
  if e, error('Wb:ERR',...
  '\ninvalid usage: single rank-2 object expected (%d)',e); end

  [q1,d1]=getQDimQS(A,1);
  [q2,d2]=getQDimQS(A,2);

  QD=[ q1, d1; q2, d2 ];
  [qd,I,d]=uniquerows(QD(:,1:end-1)); n=length(I);

  Q=qd; D=zeros(n,1);

  for i=1:n, di=QD(I{i},end); 
     if d(i)>1 && any(diff(di)), di, error('Wb:ERR',... % safeguard
     'invalid rank-2 object: severe QSpace inconsistency'); end
     D(i)=di(1);
  end

  E.Q={Q,Q};
  E.data=cell(1,n); for i=1:n, E.data{i}=eye(D(i)); end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [A,Aloc,I1,I2] = initQSpaceA0(varargin)

% init MPS block A0 representing identity (returns LRs order)
% Usage: A0= initQSpace(QL,Qs1,Qs2,...);
%
%   The tensor product of QL and Qs spans the full Hilbert space;
%   they may not be unique yet; grouping them is used to
%   determine their degeneracies.
%
% NB! returns (QL,QR=Qtot,Qs1,Qs2,...) order
%
% Wb,Sep10,06

  if nargin<2
     error('Wb:ERR','\ninvalid usage of %s\n%s', ...
     lineno('%S'), lineno('all'));
  end

  getopt('init',varargin);
     Rlast=getopt('-Rlast');
     mergeloc=getopt('mergeloc');
  varargin=getopt('get_remaining');

  nqin=length(varargin);
  for i=1:nqin
     if ~isnumeric(varargin{i}) % varargin{i}
     error('Wb:ERR','\nremaining non Q-matrix type of argument !?? %s\n%s', ...
     lineno('%S'), lineno('all')); end
  end

  QDIM=size(varargin{1},2);
  for i=2:nqin
     if size(varargin{i},2)~=QDIM
        error('Wb:ERR','\nQDIM mismatch (%d,%d)\n%s', ...
        QDIM, size(varargin{i},2), lineno('all'));
     end
  end

  if mergeloc
     for i=2:nqin
     ee{i-1}=ones(size(varargin{i},1),1); end

     for i=2:nqin
     ei=ee; ei{i-1}=varargin{i}; QQ{i-1}=mkron(ei{:}); end

     varargin={ varargin{1}, cellarrsum(QQ) };
     nqin=length(varargin);
  end

  dd=zeros(1,nqin); cn=cell(1,nqin); ee=cn; eu=cn; du=cn;
  for i=1:nqin
      dd(i)=size(varargin{i},1);
      ee{i}=ones(dd(i),1);
    % get degeneracies of input q's
      [x,x,du{i}]=uniquerows(varargin{i}); du{i}=double(du{i}(:));
      eu{i}=ones(length(du{i}),1);
  end

  QQ=cn; DD=cn;
  for i=1:nqin
    % L is fast index (just relevant within this routine)
      ei=ee; ei{i}=varargin{i}; QQ{i}=mkron(ei{:});
      ei=eu; ei{i}=du{i};       DD{i}=mkron(ei{:},'rowmajor');
    % NB! when sorting [Q1,Q2,..] below, *LAST* index will be
    % fastest through alphabetic sort -> get block dimensions
    % accordingly without `firstfast' flag (!)
  end

  DD=cat(2,DD{:}); % degeneracies of QS blocks

  [D,QDIM]=size(QQ{1}); E=speye(D,D);

  [Q1,I1,D1]=uniquerows(cat(2,QQ{:}));   m=length(D1);
  [Q2,I2,D2]=uniquerows(cellarrsum(QQ)); n=length(D2);

% determine block size Ds and DL as in D1=Ds*DL
% and mark blocks that contain non-zeros
  M=cell(m,n); mark=zeros(m,n);

  for i=1:m
     if D1(i)~=prod(DD(i,:)) error('Wb:ERR', ...
        '\nFailed to determine Q dimensions [%d; %s]\n%s', ...
         D1(i), vec2str(DD(i,:)), lineno('all'));
     end

   % NB! MatLab is col-major
   % as L is fast index -> Ls order (!)

     for j=1:n
        a=E(I1{i},I2{j}); if any(a(:))
         % NB! reshape() would not work for sparse `a' with rank>2
           mark(i,j)=1;
           M{i,j}=reshape(full(a), [ DD(i,:), D2(j)]); % LsR order
        end
     end
  end

  [I,J]=find(mark); IJ=sub2ind(size(mark),I,J);

  Q1=Q1(I,:);
  Q1=mat2cell(Q1, size(Q1,1), QDIM(ones(1,nqin)));

% need extra {} brackets to get single structure
  A=struct('Q',{{Q1{:}, Q2(J,:)}}, 'data', {M(IJ)}, 'info',{[]});

% LsR -> LRs order
  if ~Rlast
     p=1:(nqin+1); p=p([1 end 2:end-1]);
     A=permuteQS(A,p); % [1 3 2]
  end

  if nargout>1 
     if length(varargin)==2 % this includes mergeloc
          Aloc=QSpace(varargin{2},eye(size(varargin{2},1)),'operator');
     else Aloc=QSpace(varargin{2:end},'-Rlast','identity'); end
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

