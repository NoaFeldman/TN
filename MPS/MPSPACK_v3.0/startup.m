
% disp([10 'Run ~/Matlab/startup ...'])

% ------------------------------------------------------------- %
% irrelevant (matlab by default looks in current directory first!)
% path('./',path);
% r=pwd; 

  %r=getenv('MYMATLAB');
  r = '/home/noa/MPS/MPSPACK_v3.0';  % pwd;
  if isempty(r), error('Wb:ERR','MatLab path MYMATLAB not set'); end
  if ~exist(r,'dir'), r, error('Wb:ERR','invalid MatLab path MYMATLAB'); end

% adapte for working package running 'ml -mp2' 
% with path mp2 set to *this directory within cto.dat
% Wb,Oct15,14

  path(path,[r]);
  path(path,[r '/lib']);
  path(path,[r '/util']);
  path(path,[r '/bin']);
  path(path,[r '/NRG']);
  path(path,[r '/tensor']);
  path(path,[r '/setup']);
  path(path,[r '/Class']);

  clear r

% ------------------------------------------------------------- %

  ss=get(0,'ScreenSize'); ss=ss(3:4);
% fp=get(0, 'DefaultFigurePos'); fp=fp(3:4);
  fp=[590 520]; % fp=[420 400];
  fp=[ ss-fp-[4 74], fp ] ;

% starting MatLab without display sets ScreenSize = [1 1 1 1] (!)
  if ~isbatch && all(ss>1)
  set(0, 'DefaultFigurePosition',fp); end

  clear ss fp

  set(0, 'DefaultFigureName', getenv('HOST'))
  set(0, 'DefaultFigurePaperType','A4');

% works better for mjpg to avoid coarse resolution for written text
% while the rest of the figure seems good; Utopia; default: Helvetica
  set(0, 'DefaultTextFontName', 'Arial');
  set(0, 'DefaultAxesFontName', 'Arial');

  set(0, 'DefaultTextFontSize',  12 ); % default: 10
  set(0, 'DefaultAxesFontSize',  12 ); % default: 10
  set(0, 'DefaultAxesLineWidth', 1.0); % default: 0.5
  set(0, 'DefaultLineLineWidth', 1.0); % default: 0.5

% startup_aux(which(mfilename), 'opengl neverselect (see readme.txt)');
% opengl neverselect

  rand('state',sum(100*clock))

  set_global % global count (wrn,err,tst)

