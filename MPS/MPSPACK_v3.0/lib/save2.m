function save2(varargin)
% Function: save2(fname [,OPTS, usual arguments for matlab save])
%
%    Save mat file but with a few more features.
%
% Options
%
%  '-d',.. directory to save to
%  '-q'    quiet
%  '-f'    overwrite file without asking if it exists
%  '-k'    keep time stamp of existing file (if any)
%  '-x',.. list of variables to exclude
%          NB! regular expressions are supported
%          NB! *tol* would not work! use '.*tol.*' instead
%
% All remaining arguments are handed over to MatLabs save()
%
% Wb,Aug11,06  Wb,May20,07

  getopt('INIT',varargin);
     dstr =getopt('-d','');
     excl =getopt('-x',{});

         if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end

     force=getopt('-f');
     kflag=getopt('-k'); % Wb,Aug22,13

  varargin=getopt('get_remaining');

  if length(varargin)<1
     eval(['help ' mfilename]);
     if nargin, error('Wb:ERR','invalid usage'); else return; end
  end

  fname=varargin{1};
  varargin=varargin(2:end);

  if ~isempty(dstr)
     ldir=cd;
     if exist(dstr,'dir'), cd(dstr);
     else cto(dstr,'nolog'); end
  else ldir=''; end

  f=fname; m=0; nx=0; s=''; sx='';
  i=findstr(f,'./'  ); if isempty(i) && f(1)~='/', f=['./' f]; end
  i=findstr(f,'.mat'); if isempty(i), f=[f '.mat']; end

% check whether file exists (Overwrite?)
  if ~force && fexist(f), return; end

  if vflag
   % if vflag, n=char('\N'); else n='\\'; end
     wblog(1,'I/O','save data to %s',basename(f));
  end

% get full variable list via auxilliary mat-file
  if ~isempty(excl)
   % removed tmp-file usage in favor of setuser(groot,...)
   % see Archive/save2_130713.m for old version
   % Wb,Jul13,13
     evalin('caller','setuser(groot,''whos'',whos);');
     vars=getuser(groot,'whos','-rm'); ix=[];

     if ~iscell(excl), excl={excl}; end

     for i=1:length(vars)
     for q=1:length(excl)
        if ~isempty(regexp(vars(i).name,excl{q}))
         % printf('%s (%s)', vars(i).name, excl{q}); 
           ix=[ix,i]; break;
        end
     end
     end

     if ~isempty(ix) % got some variables to skip
        nx=length(ix); n=length(vars);
        sx=sprintf(' %s',vars(ix).name);
        if vflag
           if length(sx)>50, sx=[char(10) '  ' sx]; end
           fprintf('\n   skipping %d/%d variable(s): %s\n\n',nx,n,sx);
        end

        varx=vars(ix); vars(ix)=[]; m=length(vars);
        s=[' ', sprintf('%s ', vars.name)];
     end
  end

  sv=sprintf(' %s', varargin{:});
  cmd=sprintf('save %s%s%s',f,s,sv);

  if vflag>1, fprintf('\n');
     if ~isempty(varargin) || ~isempty(s)
        if length(s)>50, s=sprintf(' <%d VARS>',m); end
        fprintf(1,'   save: %s%s%s\n', f,s,sv);
     end
     fprintf(1,'   file: %s\n   path: %s\n\n', f, pwd);
  elseif vflag
     if nx>0
        if length(sx)>40, wblog(1,'I/O skipping %d vars',nx);
        else, wblog(1,'I/O skipping%s',sx); end
     end
  end

  if isbatch
     try
        evalin('caller','Imain.finished=datestr(now);');
        I=getuser(groot,'Imain'); if isfield(I,'cpuhrs')
        evalin('caller',sprintf('Imain.cpuhrs=%g;',cputime/3600-I.cpuhrs)); end
        if vflag>1, evalin('caller','structdisp(Imain);'); end
     catch l % l=lasterror; 
        fprintf(1,[...
         '\n   ERR failed to set Imain.finished time stamp' ... 
         '\n   ERR %s\n'],l.message); dispstack(l);
     end
  end

  if exist(f,'file') && kflag
     kstamp=deblank(evalc(['! date ''+%Y%m%d%H%M.%S'' -r ' f]));
  else kflag=0; end

% NB! variables larger than 2GB are not saved
% when using the standard save command % Wb,Jul13,13
% e.g. 20,000*20,000 matrix corresponds to 3.2G
% >> S.OP=randn(20000); save2 x.mat -struct S
  wsz='MATLAB:save:sizeTooBigForMATFile';
  tag=regexprep(wsz,'.*:','');

  [msg,id]=lastwarn;
  if ~isempty(findstr(id,tag)), lastwarn(''); end

  % ====================== %
    warning('off',wsz);    % NB! only turns of *display* of warning
    evalin ('caller',cmd); 
    warning('on',wsz);
  % ====================== %

  [msg,id]=lastwarn;
  if ~isempty(findstr(id,tag))
     try
      % disp(msg), disp(id) % for testing purposes // Wb,Dec20,15
        s=dbstack(1);
        if ~isempty(s)
           % s=sprintf('%s/%s:%g',s.file,s.name,s.line);
             s=sprintf(' in %s:%g',s(1).name,s(1).line);
        else s='in workspace'; end

        fprintf(1,'\n> WRN %s %s\n> WRN retry: %s -v7.3\n',...
          regexprep(msg,'cannot be saved.*','cannot be saved'),s,cmd);

        v=regexprep(msg,'.*Variable ''(.*)''.*',' $1'); fprintf(1,'\n');
        if isempty(findstr(sv,'-struct'))
           evalin('caller', ['whos ' v]);
        else
         % sv, varargin
           evalin('caller', ['whos2(''' varargin{2} ''')']);
        end; fprintf(1,'\n');

      % wblog('CMD %s\n==> retry using -v7.3 switch ...',cmd);
      % fprintf(1,'\n');

        evalin('caller', [cmd ' -v7.3']); 
     catch l % l=lasterror; 
        wblog('ERR','check save2 ...'); cmd, sv, varargin
        disp(l.message); dispstack(l);
     end
  end

  if kflag
     cmd=sprintf('touch -t %s %s',kstamp,f);
     system(cmd);
  end

  if ~isempty(ldir), cd(ldir); end

end

