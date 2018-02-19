function varargout=fixScalarOp(varargin)
% function fixScalarOp([opts,] A,B,...)
%
%   reduces "scalar" rank-3 operators to regular rank-2.
%   converting (yet keeping) the Clebsch Gordan coefficient
%   spaces to sparse.
%
% Options
%
%   '-f'  enforce reduciblity if rank-3.
%
% Wb,Aug21,13

  fflag=0; nopt=-1;
  for i=1:nargin
     if ~ischar(varargin{i}), nopt=i-1; break; end
  end
  if nopt>0
     getopt('init',varargin(1:nopt));
        fflag=getopt('-f');
     getopt('check_error');
  elseif nopt<0, error('Wb:ERR',...
    '\n   ERR invalid usage (no QSpace specified)');
  end

  mark=zeros(1,nargin);

  for i=nopt+1:nargin, n=numel(varargin{i});
     if ischar(varargin{i}), error('Wb:ERR',...
       '\n   ERR invalid usage (got intermediate option)');
     end
     for j=1:n, q=varargin{i}(j);
        if isempty(q) || numel(q.Q)<3, continue; end
        d=getDimQS(q); ok=0;
        if size(d,2)==3 && d(end,3)==1
           if isempty(q.info), q.Q(3)=[]; ok=1;
           elseif isempty(q.info.cgr), q.Q(3)=[]; ok=2;
              t=q.info.itags; if ~isempty(t)
                 if ~iscell(t), t=strread(ff,'%s','delimiter',',; '); end
                 t(3)=[]; q.info.itags=t;
              end
              varargin{i}(j)=q;
           elseif norm(q.Q{3})==0
              q.info.otype='operator'; q=QSpace(q); ok=3;
              E3=getIdentityQS(q,getvac(q));
              varargin{i}(j)=QSpace(contractQS(q,'23',E3,'12'));
            end
        end
        if ~ok && fflag
           error('Wb:ERR','\n   ERR got non-reducible QSpace (%g,%g)',i,j);
        end
     end
  end

  if nargout
     if nargin~=nargout+nopt, error('Wb:ERR',['\n   ' ... 
       'ERR invalid usage (output must match input variables, ' ...
       '%g+%g=%g)'],nargout,nopt,nargin);
     end
     varargout=varargin(nopt+1:end);
  else
     for i=find(mark), n=inputname(i);
        if isempty(n), error('Wb:ERR',['\n   ERR ' ... 
         'invalid usage (failed to access name of input variable)']); end
        assignin('caller',n,varargin{i});
     end
  end

end

