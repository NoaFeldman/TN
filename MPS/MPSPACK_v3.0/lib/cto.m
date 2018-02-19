function cto(varargin)
% Function cto - change to directory
% Usage: cto <pid>
% 
%    switches directory where pid corresponds to shortcut
%    for a directory (path ID) as specified in the cto.dat
%    
% Options
%    
%    -q   quiet
%    -p   toggle to PROJECT directory
%    
% See also unix shell script: cto <dkey>
% Wb,Jun2005  Wb,May17,06

  if nargin==0
   % eval(['help ' mfilename]);
     cd ~ ; return
  end

  getopt('init',varargin);
     vflag=~getopt({'-q','nolog'});
     pflag= getopt('-p');
  varargin=getopt('get_remaining'); narg=numel(varargin);
  
  if narg==1, pid=varargin{1}; else
     eval(['help ' mfilename]); 
     error('Wb:ERR','invalid usage'); 
  end

  if isequal(pid,'.'), return, end

  try
     if narg==1 && isequal(varargin{1},'lma')
        p=['/data/' getenv('USER') '/Matlab/Data'];
        cd(p);
     else
        error('Wb:ERR',''); 
     end

  catch l
   % hint: copy/paste through mouse seems to insert already other lines !??
     if ~isempty(find(p==10))
     wblog('ERR','directory contains NEWLINE!??'); end
     wblog('ERR','failed to change directory to >%s<\N',p);
     rethrow(l);
  end
  
  if vflag, printf('>> %s\n', pwd); end % repHome(pwd))

end

