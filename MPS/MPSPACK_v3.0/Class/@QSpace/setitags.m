function A=setitags(A,varargin)
% function A=setitags(A,{t1,t2,...} [,[i1,i2,...]])
% function A=setitags(A,i,t)
%
%    set itag(s) while preserving conj flag.
%
% adapted from QSpace/setitags.m
% Wb,Mar24,16

  if nargin<2 || nargin>3 || numel(A)~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if ~isfield(struct(A),'info') || ~isfield(A.info,'itags'), return; end

  xpat='[\*\s]';

  if nargin==3 && isnumber(varargin{1})
     i=varargin{1}; t=varargin{2};
        if ~ischar(t) || i>numel(t0)
           error('Wb:ERR','invalid usage (%g/%g)',i,numel(t0)); end
        if ~isempty(regexp(t,xpat)), error('Wb:ERR','invalid itag'); end
     w=1;
  elseif nargin>1 && iscell(varargin{1}), w=2;
     tt=varargin{1}; n=numel(tt);
        i=regexp(tt,xpat);
        if ~isempty([i{:}]), error('Wb:ERR','invalid itag(s)'); end

     if nargin>2, idx=varargin{2};
        if nargin==3 && isreal(idx) && isequal(size(tt),size(idx))
           for i=1:n
              tt{i}=sprintf('%s%g', regexprep(tt{i},'\**$',''), idx(i));
           end
        else error('Wb:ERR','\n   ERR invalid usage'); end
     end

  else error('Wb:ERR','invalid usage'); end

  for k=1:numel(A), t0=A(k).info.itags;
     switch w
        case 1
           t0{i}=set_itag(t0{i},t);

        case 2,
           if n>numel(t0),
              error('Wb:ERR','invalid itag set (%g/%g)',n,numel(t0)); end
           for i=1:n
              t0{i}=set_itag(t0{i},tt{i});
           end

         otherwise error('Wb:ERR','\n   ERR invalid switch (%g)',w);
     end
     A(k).info.itags=t0;
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function t=set_itag(t0,t)

% does not work for t0{i}='*' !)(*@ // Wb,Apr14,16
% t=regexprep(t0,'[^\*]*',t);

  t = [ regexprep(t,'\**$',''), regexprep(t0,'[^\*]*','') ];

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

