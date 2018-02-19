function mfig(varargin)
% mfig(fname [,opts])
%
%     create pdf-file with given name from current figure.
%     (default name if fname is not specified: figure.pdf)
%
% Options
%
%    'ppi',..  pixels per inch (90; larger value = smaller fig)
%              note that the actual figure size determines output figure size.
%    '-eps'    make .eps instead of .pdf
%    '-ps'     make .ps instead of .pdf
%    '-epg'    same as eps but in grayscale
%    '-png'    make .png instead of .pdf
%    '-all'    make .pdf and .eps
%    '-q'      quiet mode
%    '-f'      force overwriting of preexisting pdf-file with same name.
%
% Wb (C) 2005 / Jan26,05 / Feb06,07 / Wb,Aug02,16

  wfig={'pdf'};
  getopt ('INIT', varargin);
     ppi=getopt('ppi', 90); % 90 gets the best match in scale screen vs. figure
     force=getopt('-f');
     tflag=getopt('-t');

     for s={'-eps','-ps','-epg','-all','-png'}, s=s{1}; % Wb,Aug02,16
        if getopt(s), s=s(2:end);
           if isequal(s,'all'), wfig={'eps','pdf'}; end
           wfig={s}; break
        end
     end

     if getopt('-Q'); vflag=0;                % nothing at all
     elseif getopt({'-q','nodisp'}), vflag=1; % no display
     elseif getopt('-v'); vflag=3;            % verbose and display
     else vflag=2; % default                  % do display
     end

  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg && numel(varargin{1})==1 && isfig(varargin{1})
       fh=varargin{1}; varargin=varargin(2:end); narg=numel(varargin);
  else fh=gcf; end

  if narg
     if narg~=1 || isempty(varargin{1}) || varargin{1}(1)=='-'
        error('Wb:ERR','invalid usage'); end
     fname=varargin{1};
  else
     fname=getfigname(fh);
  end

  oo={                     % default values
    'PaperOrientation'     % 1  portrait
    'PaperPosition'        % 2  [0.25 2.5 8 6]
    'PaperPositionMode'    % 3  manual
    'PaperUnits'           % 4  inches
    'PaperSize'            % 5  [8.26772 11.6929]
    'PaperType'            % 6  A4
    'Position'             %
  };

  f0=get(fh,oo);

  fp=f0{end};
  fs=fp(3:4)/ppi; % final fig size in inch (on paper)

% [  isbatch , isequal(get(0,'ScreenSize'), [1 1 1 1]) ]
  if isbatch && isequal(get(0,'ScreenSize'),[1 1 1 1]) % && isbatch
     fs=1.5*fs;
  end

% ptype = get (fh, 'PaperType'); % save old setting (reset again below)
% set(fh,'PaperUnits','normalized');
% set(fh,'PaperType','A5');
% set(fh,'PaperPositionMode','auto') % WYSIWYG (get onscreen figure size)

% Wb,Feb24,09 - when interrupting plotting, PaperUnits may happen to be
% pixels -> way too small figure size! -> ensure inches -- tags: ppi
  set(fh,'PaperUnits','inch', ...
        'PaperOrientation','portrait',...
        'PaperSize',fs,'PaperPosition',[0 0 fs] ...
   );
  % sets 'PaperType' to '<custom>'

  h=findall(groot,'type','axes','tag','frame');;
  if ~isempty(h)
     set(h,'visible','off') % takes off frame line, yet keeps headers
  end

% wblog('TST','ppi=%g', ppi);
% fprintf(1,'\n   %18s: %s\n   %18s: %s\n   %18s: %s\n\n',...
%   'screensize',vec2str(get(0,'ScreenSize')), ...
%   'figpos', vec2str(fp), 'figsize',vec2str(fs));
% for i=1:numel(oo)
%    x={ get(fh,oo{i}), f0{i} };
%    if isnumeric(x{1})
%       x={ sprintf('[%s]',vec2str(x{1})), sprintf('[%s]',vec2str(x{2})) };
%    end
%    fprintf(1,'   %18s: %-26s %s\n',oo{i},x{:});
% end

% -------------------------------------------------------------------- %
  if vflag
     s=cell(1,2);
     if tflag, s{1}='TEST: save'; else s{1}='Save'; end
     if length(fname)<30
          s{2}=['file ' fname ' ...'];
     else s{2}=['...\n   file: ' fname]; end
     fprintf(1,['\n   %s current figure (%g,%s) to ' s{2} ...
     '\n   pwd : %s\n\n'], s{1}, get(fh,'Number'), get(fh,'Renderer'), pwd);
  end

  for i=1:numel(wfig)
      mfig_1(fh,fname,ppi,wfig{i},force,vflag,tflag);
  end
  if vflag, fprintf(1,'\n'); end

% -------------------------------------------------------------------- %
% set properties to original value

  set(fh,...
    'PaperOrientation',f0{1}, 'PaperType', f0{6},...
    'PaperUnits',f0{4}, 'PaperPosition', f0{2}   ...
  );

% NB! above resetting of paper-options ALSO CHANGES FIGURE POSITION
% when run in batch mode@!P(*Y  -- Wb,Jun12,12
  set(fh,'Position',fp);

  if ~isempty(h), set(h,'visible','on'); end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function mfig_1(fh,fname,ppi,fext,force,vflag,tflag)

  eps=0;
  switch fext
    case 'eps', ffmt='epsc2';           eps=1;
    case 'epg', ffmt='eps'; fext='eps'; eps=1;
    case  'ps', ffmt='epsc2';           eps=1;
    case 'pdf', ffmt=fext;
    case 'png', ffmt=fext;
    otherwise fext, error('Wb:ERR','\n   ERR invalid fext');
  end

  if isempty(regexpi(fname,['.' fext '$'])), fname = [ fname '.' fext ]; end
  if isempty(fileparts(fname)), fname=['./' fname]; end

% wblog(' * ','%s %s %g %g',fext,ffmt,vflag,tflag); 

  if ~force 
     [i,fname]=fexist(fname,tflag);
     if i, % disp(' huh? ')
     return, end
  end

  set(fh,'FileName',fname);

  if vflag>1 || tflag
     fprintf(1,'   saveas(fh,''%s'',''%s'');\n',fname,ffmt);
   % fprintf(1,'   %s\n',CMD);
     if tflag, return; end
  end

% wblog('TST','saveas(%g,%s,%s) ...',fh,fname,ffmt);
% THIS CAN TERMINATE MATLAB INSTANTANEOUSLY WITHOUT ANY FURTHER MESSAGE
% use direct command, instead [eval: matlab-bug] // Wb,Feb18,11
% CMD=sprintf('saveas(fh, ''%s'', ''%s'');', fname, ffmt);
% eval(CMD);
  saveas(fh,fname,ffmt);
% wblog('TST','ok');

  if vflag>1
     if eps, eval(['! gv --watch ' fname ' &']);
     else
      % p='kpdf'; % acroread
      % eval(['! ' p ' ' fname ' &']);
        eval(['! vi.pl ' fname ' &']);
     end
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

