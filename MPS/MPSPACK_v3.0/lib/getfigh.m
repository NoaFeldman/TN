function h=getfigh(tag,varargin)
% Function h=getfigh(tag,figfile [,OPTS])
%
%    Get axis set with given tag
%    either from current Matlab environment
%    or if it does not exist there from given figfile.
%
% Wb,Feb21,08

  getopt('init',varargin);
     ezflag=getopt('-ez');
  varargin=getopt('get_remaining');

  narg=length(varargin);
  if narg
     if narg~=1 || ~ischar(varargin{1})
     error('Wb:ERR','invalid usage'); end
     figfile=varargin{1};
  else figfile=''; end

  if narg
   % set figure.tag to filename for subsequent calls to allow to 
   % check whether figure is already open
     ftag=['gfh_' figfile];
     f=findall(groot,'type','figure','tag',ftag);
     if isempty(f)
        wblog('I/O','loading %s ...',figfile);
        f=open(figfile); set(f,'Visible','off','tag',ftag);
     end

     if isempty(f), error('Wb:ERR','invalid usage (figure not found)'); 
     elseif numel(f)>1
        wblog('WRN','more than one figure set found (%g)',numel(f))
        f=f(1);
     end

     h=findall(f,'type','axes','tag',tag);
  else
     h=findall(groot,'type','axes','tag',tag);
  end

  if isempty(h)
     s=sprintf('failed to find / open axes %s (%s)',tag,figfile);
     if ezflag, wblog('Wb:ERR',s);
     else error('Wb:ERR',s); end
  end

end

