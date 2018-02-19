function P=PSet(varargin)
% Class: P=PSet([P,] 'B',bb, 'T',tt, ...)
%
%   automated parameter set handling for multiple parameter set.
%   If first argument is already a PSet, update it accordingly.
%
%   NB! parameters are traversed col-major
%   (first index / parameter is fastest).
%
% Wb,Jul31,08

  pflag=0;
  if nargin, P=varargin{1};
     if isa(P,'PSet')
      % may be used to update / add to input PSet
        P=struct(P); pflag=1;
        varargin=varargin(2:end);
     elseif builtin('isfield',P,'name') && builtin('isfield',P,'data')
      % initialization from P structure
        P=struct('name',{P.name},'data',{P.data});
        varargin=varargin(2:end);
     else
        P=struct('name',{{}},'data',{{}});
     end
  else
   % construct from scratch
     P=struct('name',{{}},'data',{{}});
  end

  narg=length(varargin);

  if mod(narg,2)==1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end
  for i=1:2:narg
      if ~ischar(varargin{i}) || ~isnumeric(varargin{i+1})
      error('Wb:ERR','invalid PSet constructor'); end
  end

  P.name={P.name{:}, varargin{1:2:end}};
  P.data={P.data{:}, varargin{2:2:end}};

  [x,I,d]=uniquerows(strvcat(P.name{:}));
  if any(d>1)
     for i=1:length(I)
         if length(I{i})>1, wblog('WRN',...
         'overwriting existing parameter set ''%s''',P.name{I{i}(1)}); end

         I{i}=I{i}(1:end-1);
     end
     I=cat(2,I{:}); P.name(I)=[]; P.data(I)=[];
  end

  r=length(P.data); s=zeros(1,r);
  for i=1:r, s(i)=prod(size(P.data{i})); end

  P.r=r;
  P.n=prod(s);
  P.s=s;
  P.i=0;
  P.t=[];

  if isempty(s), P.n=0; end

  P=class(P,'PSet');

  if pflag && ~nargout, n=inputname(1);
     if ~isempty(n), assignin('caller',n,P); clear P; end
  end

end

