function rval = getopt(varargin)
% function rval = getopt(varargin)
%
%    options can be grouped (i.e. multiple options implying the same)
%    getopt({opt1,opt2} [,default_value])
%
% Wb,Nov05  Wb,Jan09,08

  persistent args errcount casesensitive

  if nargin<1
  eval(['help ' mfilename]); return; end

  if nargout, rval=[]; end

% wblog('nargin = %d (%s)', nargi, varargin{1}); args

  if isnumeric(varargin{1})
     lflag=varargin{1}; % log flag
     varargin=varargin(2:end);
  else lflag=0; end

  nargi=length(varargin);

  if ~iscell(args), args = {args}; end

% INIT ----------------------------------------------------------
  if isequal(varargin{1}, 'init')
     args = varargin{2:end};
     errcount=0;
     casesensitive=0;

  elseif isequal(varargin{1}, 'INIT')
     args = varargin{2:end};
     errcount=0;
     casesensitive=1;

% CHECK IF SOMETHING IS LEFT ------------------------------------
  elseif strcmpi('check_error', varargin{1}) & nargi==1
     if length(args)>0
      % wblog('WRN - the following options will be discarded ...');
        wblog(2,'ERR','\Ninvalid options ...');
        args  % args{1}
        errcount = errcount + 1;

        fprintf(1,'Stack: '); s='-> ';

        [stack, index] = dbstack; m=length(stack);
        for i=m:-1:1
           if i==1, s=''; end
           fprintf(1,'%s::%d %s', stack(i).name, stack(i).line, s);
        end

        if length(stack), fprintf(1,'\n\n\n'); end

        error(' ');
      % error('wb:getopt','invalid options');
     end

     if nargout, rval=errcount; end
     args={}; errcount=0;

% CHECK IF SOMETHING IS LEFT ------------------------------------
  elseif strcmpi('get_remaining', varargin{1}) & nargi==1

     rval=args;
     errcount=0;

  elseif strcmpi('get_last', varargin{1}) & nargi<=2

     n=length(args);
     if n==0
        if nargi>1, rval=varargin{2};
        else error('Wb:ERR','invalid usage (missing default value)');  end
     else
        if n==1, rval=args{1};
        else getopt('check_error'); end
     end

% SEARCH FOR VALUE ----------------------------------------------
  elseif nargi<=2    %% vname [,default]
     NAME=varargin{1}; narg=length(args); found=0;

     if ~iscell(NAME), NAME={NAME}; end
     no=length(NAME);

     for io=1:no
        if casesensitive, vname=NAME{io}; else vname=lower(NAME{io}); end
        for i=1:narg, if ~ischar(args{i}), continue; end
            if casesensitive, ai=args{i}; else ai=lower(args{i}); end
            if isequal(vname,ai), found=1; break; end
        end
        if found, break, end
     end

     if nargi<2
        if found, if lflag, wblog(1,' * ','%s',NAME{io}); end
             rval=1; args(i)=[];
           % accept flag with subsequent number
           % e.g. may be last argument, to be extracted using 'get_last'
           % Wb,Apr24,12
             if i<=length(args) && isnumeric(args{i}) && NAME{io}(1)~='-'
                if isscalar(args{i}), wblog('WRN',...
                  'asking for optional flag, yet value %g specified !??',args{i});
                   rval=(args{i}~=0); args(i)=[];
                else error('Wb:ERR',...
                'asking for optional flag, yet array specified !??'); end
             end
        else rval=0; end
        return
     else
        if found
           if i==narg
              wblog('ERR - value expected for option %s', vname);
              errcount=errcount+1;
              args(narg)=[];
           else
              rval=args{i+1}; args(i:i+1)=[];
              if lflag
                 if isnumeric(rval)
                      wblog(1,' * ','%-8s: %g',vname,rval);
                 else wblog(1,' * ','%-8s: %s',vname,rval); end
              end
           end
        else rval= varargin{2}; end
        return
     end
  else
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

% if nargout==0; clear rval; end

end

