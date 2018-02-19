function lma(varargin)
% Function: lma - load .mat file from local directory (LMU net)
% Usage: lma ([OPTS, <file pattern>])
%
%    <fpat>  file pattern [*.mat]
%
% Options
%
%   -a       show all mat files
%   -n, ..   number of records on display
%   -d,'..'  specifies directory tag as recognized by cto (default directory: lma)
%   '.'      choose local directory (same as -d .)
%
% Examples
%
%    lma -k 12 . kro*    % search current directory for at max 12 kro*.mat files
%
% Wb,Oct25,05  Wb,May11,07

  mdir='lma'; pat='*.mat';

  getopt('init',varargin);
     k    =getopt('-n',4);
     aflag=getopt('-a');  if aflag, k=Inf; end
     mdir =getopt('-d','lma');
     ldir =getopt('.');   % choose local directory
  varargin=getopt('get_remaining');

  cdir=cd;
  try if ~ldir, cto(mdir); end

    if length(varargin)
       pat=varargin{1}; varargin=varargin(2:end);
       if all(isalpha(pat)), pat=[ '*' pat '*.mat' ]; end
    end

    d=dir2(pat,varargin{:}); n=length(d);
    if isempty(d)
       fprintf(1,'   No files "%s" found.\n',pat);
       cd(cdir); return
    end

    [it,ii]=sort(datenum(cat(1,d(:).date),'yyyy-mm-dd HH:MM:SS')); d=d(ii);
    ii=max([1,n-k]):n;

    fprintf(1,'\n     #   DATE                    SIZE     FILENAME\n');
    for i=ii
        s=num2str2(d(i).bytes,'-bytes');
        fprintf(1,'    %2d   %s  %10s  %s\n', ...
        ii(end)-i+1, d(i).date, s, basename(d(i).name) );
    end; fprintf(1,'\n');

    if length(it)>1 
       i=input(sprintf('>> Select file <#[1]> (0=none) ... '));
       if isempty(i), i=1; end % i=ii(end);
    else
     % i=input(sprintf('\n>> Select file <#|[0=none]> ... '));
     % if isempty(i), i=0; end % i=ii(end);
       i=1;
    end

    if i<1
     % fprintf(1,'>> None.\n\n')
       fprintf(1,'\n'); cd(cdir); return
    end
    i=ii(end-i+1);

  % fprintf(1,'\nSelected: %s', d(i).name);

    fprintf(1,'>> load %s\n', d(i).name);
    evalin('caller', ['load(''' pwd '/' d(i).name ''')']);

  catch l % l=lasterror; 

    cd(cdir); disp(l.stack(1));
    if length(l.stack)>1, disp(l.stack(end)); end
    rethrow(l)

  end

  cd(cdir);

end

