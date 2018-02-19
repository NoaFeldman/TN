function varargout=fixAbelian(varargin)
% function fixAbelian(A,B,...)
% Options
%
%    -xop3    reduce "scalar" rank-3 operators to regular rank-2.
%
%    -addCGS  add (trivial) CGS space to all-abelian operators.
%    -op3     together with -addCGS: make operator of rank 3 if rank 2.
%
% Wb,Feb20,13

  isopt=zeros(1,nargin);
  for ia=1:nargin
     if ischar(varargin{ia}), isopt(ia)=1; end
  end

  i=find(isopt); IA=1:nargin; addcgs=0;
  if ~isempty(i)
     opts=varargin(i); IA(i)=[]; nargs=numel(IA);
     getopt('init',opts);
        xop3=getopt('-xop3'); if ~xop3
        addcgs=getopt('-addCGS');
        makeop=getopt('-op3'); end
     getopt('check_error');
  else xop3=0; nargs=nargin; end

  if nargout
     if nargs~=nargout, error('Wb:ERR',...
       '\n   ERR invalid usage (must match input to output variables)'); end
  end

  mark=zeros(size(IA));

  if addcgs==0
     for ia=1:nargs
        A=varargin{IA(ia)}; if ~isempty(A), isa=isAbelian(A); isa=isa(:);
        if any(isa)
           for k=1:numel(A)
              if isa(k) && xop3
                 d=getDimQS(A(k));
                 if length(d)==3 && d(3)==1, A(k).Q(3)=[]; end
              end
              if isa(k)>1, A(k).info={}; end
           end
           varargin{IA(ia)}=A; mark(ia)=1;
        end, end
     end
  else
     I=struct('qtype','','otype','','cgs',{{}},'itags',{cell(1,0)});
     for ia=1:nargs
        A=varargin{IA(ia)}; m=0;
        if isempty(A) || ~all(reshape(isAbelian(A),[],1))
           error('Wb:ERR','\n   ERR got non-abelian symmetries'); 
        end

        for k=1:numel(A)
           if isempty(A(k).info), A(k).info=I; m=1; end
           if makeop && numel(A(k).Q)==2, m=1;
              A(k).Q{3}=A(k).Q{1}-A(k).Q{2};
              A(k).info.otype='operator';
              A(k).info.cgs=[];
           end

           s=size(A(k).Q{1}); r=numel(A(k).Q);

           if isempty(A(k).info.qtype), m=1;
              q=repmat(',A',1,s(2));
              A(k).info.qtype=q(2:end);
           end

           if isempty(A(k).info.cgs), m=1;
              if r==2
                   A(k).info.cgs=repmat({1},s);
              else A(k).info.cgs=repmat(struct('S',ones(1,r),'data',1),s);
              end
           end
        end
        varargin{IA(ia)}=A; mark(ia)=m;
     end
  end

  
  if nargout, varargout=varargin(IA);
  else
     for ia=find(mark), k=IA(ia);
        n=inputname(k); if isempty(n), error('Wb:ERR',['\n   ERR ' ... 
         'invalid usage (failed to access name of input variable)']); end
        assignin('caller',n,varargin{k});
     end
  end

end

