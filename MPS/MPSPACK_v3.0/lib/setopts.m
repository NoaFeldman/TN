function oo=setopts(varargin)
% Function: setopts(olist, opt1, 'opt2', ...) 
%
%    Adds options to given option list olist; if no output
%    variable is specified and olist is specified on input,
%    output directly replaces olist.
%
% Adding options with values
%
%    opt1, ...     adds opt1 if set in calling routine
%   'opt1:', val   adds opt1 with value specified
%   {'opt1',val}   adds opt1 with value specified
%   'var1?'        adds option var1 with value of variable var1 in caller space
%   {'opt1'}       removes options opt1 from olist if it exists
%
% Adding flags without values
%
%    '--flag?'      adds flag '-flag' to olist iff variable flag exists and is ~=0
%    '--flag'       adds flag '-flag' to olist 
%    '~-flag'       removes flag '-flag' from olist
%
% Wb,Oct23,06  Wb,Jun14,07

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % NB! See also MatLab's join() routine! %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  evalin('caller','global gopt_val__ gopt_flag__');
  global gopt_val__ gopt_flag__

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

% NB! order of input arguments is crucial! (using inputname below)
% ==> qflag may only appear at the very first position
  QFLAG=0; narg=length(varargin);
  
  vn=cell(1,narg);
  for i=1:narg, vn{i}=inputname(i); end
  
% check for -q flag (accept anywhere)
  for i=1:narg, q=varargin{i};
     if ischar(q) && ~isempty(regexp(q,'^-[qQ]$'))
        switch q
           case '-q', QFLAG=2;
           case '-Q', QFLAG=2;
           otherwise, error('Wb:ERR','invalid usage');
        end
        varargin(i)=[]; vn(i)=[]; break;
     end
  end

  oo=varargin{1}; varargin=varargin(2:end);
  on=vn{1}; vn=vn(2:end);

  if isequal(oo,'-'), oo={}; noassign=1; else noassign=0; end

  n=length(varargin); i=0;
  while i<n % keep flexible: flags vs. options with values

     i=i+1; nm=varargin{i}; gotval=0; qflag=QFLAG;
     if iscell(nm) && isempty(vn{i}), l=numel(nm);
      % value is specified explicitely by {'opt,val}

        if l==2,     val=nm{2}; nm=nm{1}; gotval=1;
        elseif l==1, val=[];    nm=nm{1}; gotval=-1; % remove from options!
        else error('Wb:ERR','invalid usage {''name'',val }'); end

     elseif i<n && ~iscell(varargin{i+1}) && ~ischar(varargin{i+1}) && ...
        isempty(vn{i+1})
        i=i+1; val=varargin{i}; gotval=1;

     else %if ~ischar(nm) % variable may carry value of type string!

      % option name = name of variable on input
        if ~isempty(vn{i}), val=nm; nm=vn{i}; gotval=1; end
      % else gotval=0.
     end

     if ~ischar(nm) || isempty(nm) error('Wb:ERR',...
     'invalid usage (option name not of type char)'); end
      
   % alternative way specify value explicitely: 'opt:',val, ...
     if ~gotval
        if nm(end)==':', i=i+1;
           if i>n, error('Wb:ERR','missing value for `%s'')',nm); end
           nm=nm(1:end-1); val=varargin{i}; gotval=1;
        elseif nm(1)=='-' || nm(1)=='~' % got option/flag without value
           gotval=nm(1); nm=nm(2:end);
           if nm(end)=='?', cflag=1; nm=nm(1:end-1); else cflag=0; end
        end
     end

     if ~gotval % get value from workspace
        switch nm(end)
           case '?', qflag=1; nm=nm(1:end-1);
           case '!', qflag=2; nm=nm(1:end-1);
           otherwise qflag=0;
        end

        evalin('caller',sprintf([
          'gopt_flag__=0; if exist(''%s'',''var'')==1, ' ...
          'gopt_val__=%s; gopt_flag__=1; end'],nm,nm));
        if gopt_flag__, gotval=1; val=gopt_val__; end

        if qflag<QFLAG, qflag=QFLAG; end
        if qflag && ~gotval, continue; end
     end

     k=findcell(oo,nm);

     if length(k)>1, wblog('WRN',...
     'option `%s'' exists %d times !?',nm,length(k)); end

   % check option/flag without value
     if gotval=='-', gopt_flag__=1;
        if cflag, s=nm; s=regexprep(s,'^[-]+','');
           cmd=sprintf(['if ~exist(''%s'',''var'') || ' ...
             'numel(%s)~=1 || ~%s, gopt_flag__=0; end'],s,s,s);
           evalin('caller',cmd);
        end
        if isempty(k) && gopt_flag__, oo{end+1}=nm; end
        continue
     elseif gotval=='~'  % unset option without value
        if ~isempty(k)
           if ~qflag, wblog(' * ','removing option `%s''',nm); end
           oo(k)=[];
        end
        continue
     end

     if gotval<=0, if gotval==0
        error('Wb:ERR','missing value for option ''%s''',nm); end

        if ~isempty(k), oo([k,k+1])=[]; k=[]; end
        continue
     elseif length(k)>1
        oo([k(2:end), k(2:end)+1])=[]; k=k(1);
     end

     if ~isempty(k) % option already exists
        if qflag~=2 && ~isequal(oo{k+1},val) % if gotval
           wblog('<o>','Replace existing option ''%s''',nm); % end

           if (ischar(oo{k+1}) || numel(oo{k+1})<12) && ...
              (ischar(val) || numel(val)<12)
              s=evalc('disp(oo{k+1})');
              fprintf(1,'  %s\n->%s', s(1:end-2), evalc('disp(val)'))
           end
        end
     else
        k=length(oo)+1;
     end

     oo(k:k+1)={nm,val};
  end

  if nargout==0 && ~noassign
     if isempty(on), error('Wb:ERR','invalid usage'); end
     assignin('caller',on,oo); clear oo
  end

  clear global gopt_val__ gopt_flag__

end

