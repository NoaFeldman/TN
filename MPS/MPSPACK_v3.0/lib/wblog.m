function wblog(varargin)
% wblog - logging routine
%
%    Usage : wblog(fmt, ...)
%    Output: filename:linenr (output as with sprintf())
%
% format string
%
%    \n  leads to newline with usual file:line header
%    \N  leads to old fashioned newline, i.e. with no header
%
%    except for the very first appearances of \N or \n in fmt,
%    which are interpreted as global leading \n
%
% Wb,Nov18,05: allow to redirect output to file
%
%    wblog('setfid',fid)     - redirects the following wblog's to file fid
%    wblog('setfout',fname)  - redirects the following wblog's to file fname
%    wblog('unsetfid')
%    wblog('unsetfout')      - unset the output to file
%
% Based on lineno() from C.Denham which uses dbstack()
% for more info: see help to lineno()
%
% Wb,Oct08,03

  persistent fid newfile llday

  if nargin==0, eval(['help ' mfilename]); return; end
  if nargin==1

     if isequal(varargin{1},'-ping')
        llday=check_nextday(llday,fid);
        return
     end

     if isequal(varargin{1},'unsetfid')
     fid=[]; return, end

     if isequal(varargin{1},'unsetfout')
        if ~newfile
           wblog('ERR','fid set by calling routine.');
        elseif ~isempty(fid)
           fclose(fid);
        else wblog('ERR','No file to close.'); end

        fid=[]; return
     end

  elseif nargin==2

     if isequal(varargin{1},'setfid')

        if ~isempty(fid)
        wblog('ERR','must close active fid first.'); return; end

        fid=varargin{2}; newfile=0;
        try, fprintf(fid,'');
        catch,fid=[]; wblog('ERR','Invalid fid.'); end
        return

     elseif isequal(varargin{1},'setfout')

        if ~isempty(fid)
        wblog('ERR','must close active fid first.'); return; end

        fname=varargin{2};
        fid=fopen(fname,'w'); newfile=1;
        try, fprintf(fid,'');
        catch,fid=[]; wblog('ERR','Cannot write to file %s.',fname); end
        return

     end
  end

  if ~ischar(varargin{1}) & length(varargin{1})==1
       idx=varargin{1}; varargin={varargin{2:end}};
  else idx=0; end

  s=lineno_aux(2+idx);
  t=datestr(now);

  tag=''; TAG=tag;
  if length(varargin)>1 && isempty(findstr(varargin{1},'%'))
     tag=varargin{1}; varargin={varargin{2:end}};
     TAG=tag; cstflag=0;
     switch tag
       case {'LL1'}, tag='.  ';
       case {'LL2'}, tag=' . ';
       case {'LL3'}, tag='  .';
       case {'CST'}, cstflag=1; tag='';
       case {'Error'}, tag='ERR';
       case {'Warning'}, tag='WRN';
     end
     if ~cstflag, tag=sprintf('%3s ',tag); end

     if ~isempty(tag)
        u=get(0,'UserData'); f=deblank(lower(tag));
        if isempty(u) || ~isfield(u,'err_count')
        set_global; u=get(0,'UserData'); end

        g=u.err_count;
        if isfield(g,f)
           g=setfield(g,f, getfield(g,f)+1);
           u.err_count=g;
           set(0,'UserData',u);
        end
     end
  end

  hstr=sprintf('%-20.20s %s  %s', shortfstr(s,20), t(13:end), tag);

  fstr=hstr; % fstr(:)=' '; % fstr(end-1)='.'; % follow up string (newlines)

% substitute line breaks
  fmt = varargin{1};

  for i=1:length(fmt), if fmt(i)~=10, break, end, end
  if i>1, inl(i-1); fmt=fmt(i:end); end

  while length(fmt)>2 && (isequal(fmt(1:2),'\N') || isequal(fmt(1:2),'\n'))
     fprintf(1,'\n'); fmt=fmt(3:end);
  end

  fmt = strrep(fmt,char(10),[10 fstr]);
  fmt = strrep(fmt, '\n', [10 fstr]);
  fmt = strrep(fmt, '\N', char(10));

  if     length(fmt)>1 && isequal(fmt(1:2),'\\'), fmt=fmt(3:end);
  elseif length(fmt)>1 && isequal(fmt(1:2),'\r'), fmt=['\r' hstr fmt(3:end)];
  else fmt=[hstr fmt]; end

% tags: no final CR \n \b \r
  if isequal(fmt(max(1,end-1):end),'\\'), fmt=fmt(1:end-2);
  else fmt=[fmt '\n']; end

% log day switch
  llday=check_nextday(llday,fid);

% print the whole thing
  if isempty(fid)
     fprintf(1,fmt,varargin{2:end}); % 1=stdout
     if length(TAG)>3
        switch TAG % msg identifier must have structure 's1:s2:s3'
          case 'Warning', inl(1); warning('Wb:lib:WRN',''); inl(1);
          case 'Error',   inl(2); error  ('Wb:lib:ERR','');
        end
     end
  else
     fprintf(fid,fmt,varargin{2:end});
  end

  if ~isempty(who('global','DEBUG'))
      global DEBUG
      if DEBUG, switch TAG
       case 'WRN', inl(1); dbstack; beep; inl(1); keyboard
       case 'ERR', inl(2); dbstack; beep; inl(1); keyboard
      end, end
  end

end

% -------------------------------------------------------------------- %
% NB! since wblog() is called within getopt avoid
% wblog() to use routines with varargin calling getopt()
% -------------------------------------------------------------------- %

function line = lineno_aux(id)
% short version of lineno()

  [stack,index]=dbstack;
  stack=stack(min([1+id, length(stack)]));

% default format: file>name:line
% but drop name if equal to file (i.e. the actual function)
  if ~isequal(stack.file, [stack.name '.m']);
       line = sprintf('%s>%s:%d', stack.file, stack.name, stack.line);
  else line = sprintf('%s:%d', stack.file, stack.line); end
% line=shortfstr(line,len);

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% check day of last log

function llday=check_nextday(llday,fid)

% d=str2num(datestr(now,'dd')); % check only day (2 digits)
  d=datevec(now); d=d(3); % [year,month,day,hour,min,sec]
  if ~isequal(d,llday)
     if ~isempty(llday)
        if isempty(fid), f=1; else f=fid; end
        fprintf(f,'\n\n>> TODAY %s :: (%d->%d)',datestr(now),llday,d);
      % kinit -v 
        fprintf(f,'\n\n');
     end
     llday=d;
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
