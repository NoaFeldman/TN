% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% ==================================================================== %
% -------------------------------------------------------------------- %

  getopt('init',varargin);
     vflag=getopt('-v');
  getopt('check_error');
  varargin=getopt('get_remaining'); narg=length(varargin);
  a=getopt('get_last',[]);

% fprintf(1,'\n'); eval(['help ' mfilename]);
  if length(varargin)<2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if nargin>1 && isnumeric(varargin{1})
     k=varargin{1}; varargin=varargin(2:end);
     narg=length(varargin);
  end

% eval(['help ' mfilename]);
  if nargin==0, helpthis; return; end

  P=PSet('epsd',...);
  while next(P), [p,pstr]=P(); structexp(p);
     ...
  end

  for ip=1:P.n, [p,pstr,tstr]=P(ip); structexp(p);
     if P.n>1, banner('%s\n%s',pstr,tstr); end
     ...
  end

  [it,nt]=P(ip,'T');

  drawnow; k=keyiter(k,n); if k<1 || k>n, break; end
  k=[1 15]; while keyiter(k), ...; end

  if ~nargout
     n=inputname(1); if isempty(n)
     error('Wb:ERR','failed to obtain inputname of input argument'); end

     assignin('caller',n,S);
     clear S
  end


ah=smaxis(2,2,'tag',mfilename);
setax(ah(1,1))

  varargout=cell(1,max(1,nargout));
  [varargout{:}]=hasQOverlap(QSpace(varargin{1}), varargin{2:end});

  try
    djdjddj
  catch l
     wblog('gotcha'); % l=lasterror; 
     dispstack(l); rethrow(l);
  end

  datestr(now,'yymmdd')

% tags: color
  set(gcf,'DefaultAxesColorOrder', [
     .0  .0  .0
      1  .5   0    % orange
     .1  .5   1    % blue
     .88 .68 .34   % orange.in
     .13 .61 .28   % green.in
     .27 .34 .44   % darkblue.in
      1  1 .7      % marker yellow
  ]);


% -------------------------------------------------------------------- %
% interactive

  i=input('question <[1]|0> '); % shown without extra newline

  if isempty(i); i=1;
  elseif i~=1; return; end

  if ~isbatch
     q=input([
      '\n   WRN sure to (re)run fdm/NRG?' ...
      '\n   WRN press return to continue / Ctrl-C to terminate ... \n\n'
     ],'s');
  end

% -------------------------------------------------------------------- %

