function rload(varargin)
% Usage: rload(['-f',] rhost,fname)
% 
%    option -f forces to copy the file
%    (even if it (already) exists locally)
% 
%    if no path is specified with fname
%    default path on rhost is assumed to be $LMA
% 
%    file will be copied to localhost//$LMA
% 
% See also loadfs()
% Wb,Oct29,06

  fflag=0;

  if nargin
     if isequal(varargin{1},'-f'), fflag=1; varargin=varargin(2:end); end
  end

  if length(varargin)~=2, eval(['help ' mfilename]); return; end

  rhost=varargin{1};
  fname=varargin{2};

  lma='/data/$USER/Matlab/Data/';

  if length(fname)>4 && ~isequal(fname(end-3:end),'.mat')
  fname=[fname '.mat']; end

  ii=findstr(fname,'/');
  if isempty(ii)
     floc=[lma rhost(1:3) '_' fname];
     fname=[lma fname];
  else
     floc=[lma rhost(1:3) '_' fname(ii(end)+1:end) ];
  end

  if isequal(hostname,rhost)
     printf('\n   host is local (%s)', rhost);
     printf(  '   loading file %s ...\n', strrep(fname,lma,''));
     evalin('caller', ['load(''' fname ''')']);
     return
  end

  cmd{1} = ['!scp -p $USER@' rhost ':' fname ' ' floc ];
  cmd{2} = ['load(''' floc ''')'];
% cmd{2} = ['load(''' floc '''); !rm ' floc];

  try
     if fflag || ~exist(floc,'file')
        printf('\n   copying file to LMA ...\n');
        eval(cmd{1});
        printf('\n   loading file %s ...\n', strrep(floc,lma,''));
     else
        printf('\n   file exists; taking local copy.');
        printf(  '   loading file %s ...\n', strrep(floc,lma,''));
     end

     evalin('caller', cmd{2}); 
  catch l
     disp(cmd{1})
     disp(cmd{2})
     rethrow(l)
  end

end


