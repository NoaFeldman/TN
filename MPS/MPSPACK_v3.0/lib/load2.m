function S=load2(varargin)
% Function: S=load2([OPTS,] fname [,vars])
%
%   Similar to load, but with extra features.
%
%   load2('var1?',var2,var3)  % load var1 only if it exists in file
%   S=load2(..)               % return loaded data in data structure
%
% Options
%
%  '-q'     quiet (do not perform warning if variables were not found)
%  '-d',..  directory where mat file is to be loaded from
%           accepts shortcuts as in cto.
%  '-x',..  exclude variables with given pattern (also applies to Class name)
%  'smax',. max. variable size to include (ignores 'var?' type input)
%
% Wb,Nov02,05  Wb,May24,05

  getopt('init',varargin);
     dstr =getopt('-d',[]);
     vflag=getopt('-v'   );
     qflag=getopt('-q'   );
     xpat =getopt('-x',[]);
     smax =getopt('smax',Inf);
  varargin=getopt('get_remaining');

  if length(varargin)<1
     eval(['help ' mfilename]);
     if nargin, error('Wb:ERR','invalid usage'); else return; end
  end

  fname=varargin{1}; varargin=varargin(2:end);

  if ~ischar(fname), error('Wb:ERR','invalid file name'); end

  if ~isempty(dstr), ldir=pwd;
     if exist(dstr,'dir'), cd(dstr);
     else cto(dstr,'nolog'); end
     if isempty(findstr(fname,'./')), fname=['./' fname]; end
  end

  if ~exist(fname,'file') && ~isempty(findstr(fname,'*'))
     q=dir(fname); nq=numel(q); fname_=fname;
     if nq==1
        fname=q.name;
     elseif nq>1, disp({q.name}')
        error('Wb:ERR','\n   ERR multiple files found'); 
     end
  end

  if ~exist(fname,'file')
     f2=[fname '.mat'];
     if ~exist(f2,'file')
         wblog('ERR','file not found\N\N   file: %s\N   pwd : %s\N', fname, pwd);
         error('Wb:ERR','file not found');
     end
     fname=f2;
  end

% show size if file is larger than 9MB
% NB! isempty(fs) may happen if fname is not in local directory
  fs=dir(fname); if isempty(fs), fs=dir(which(fname)); end; s=fs.bytes;
  if fs.bytes<9E6 && ~vflag || qflag, s='';
  else
     s=sprintf(' [%s, %s]',datestr(fs.datenum,'dd/mm/yyyy'), ...
     num2str2(fs.bytes,'-bytes','fmt','%.0f'));
  end

  if qflag
     warning off MATLAB:load:variableNotFound;
     if ~isempty(find(cat(2,varargin{:})=='?'))
     disp(strhcat(varargin))
     wblog('WRN','load2 with option -q must not have var? arguments');

     end
  elseif ~isempty(xpat) || smax<Inf
     S=whos('-file',fname);
     m=length(S); mark=zeros(1,m);
     for i=1:m
        if S(i).bytes>smax || ~isempty(xpat) && (...
           ~isempty(regexp(S(i).name, xpat)) || ...
           ~isempty(regexp(S(i).class,xpat)) )
           mark(i)=1;
        end
     end
     I=find(mark); if ~isempty(I)
        if ~qflag
        wblog('<i>','skipping %g variables\n%s',...
          length(I), sprintf('%s ',S(I).name)); end
        S(I)=[];
     end
     varargin={S.name};
  else
   % check optional variables of type 'var?'
     m=length(varargin); mark=zeros(1,m);
     for i=1:m
        if isequal(varargin{i}(end),'?')
           varargin{i}=varargin{i}(1:end-1);
           mark(i)=1;
        end
     end

     I=find(mark); if ~isempty(I)
        S=whos('-file',fname,varargin{I});
        getopt('init',{S.name});
           for j=1:length(I)
           if getopt(varargin{I(j)}), mark(I(j))=0; end, end
        getopt('check_error');
        I=find(mark); if ~isempty(I)
         % wblog('TST','vars { %s } not found.',strhcat(varargin(I)));
           varargin(I)=[];
        end
     end
  end

  if ~nargout && ~qflag || vflag % || length(varargin)==0
     n=which(fname); if isempty(n), n=which([fname '.mat']); end
     if isempty(n), n=fname; end; n=repHome(n);
     if length(n)>30, [d,f,x]=fileparts(n);
        if length(d)<10
         % wblog(1,'I/O','loading file from %s%s',d,s); s='';
           fprintf(1,'\n   loading file from %s%s',d,s); s='';
        else
           wblog(1,'I/O','loading file ...');
           fprintf(1,'\n   dir: %s',d);
        end
        fprintf(1,'\n   mat: %s%s\n\n',[f x],s);
     else
        wblog(1,'I/O','loading file%s %s',s,n);
     end
  end

  if ~nargout
     evalin('caller', ['load ' fname sprintf(' %s',varargin{:})]); 
  else
   % return structure (only when variables are named explicitely)
     S=load(fname,varargin{:});
  end

  if qflag
  warning on MATLAB:load:variableNotFound; end

  if ~isempty(dstr), cd(ldir); end

  if ~nargout, clear S; end

end

