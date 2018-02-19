function meps(fname)
% meps([fname='figure.eps'])
%
%     make (e)ps file from current plot window
%     default output file: figure.eps
%
% Wb (C) 2001  (Apr22,01)

% when called like 'meps fname' then all the input args are strings!

  if ~nargin, fname=getfigname(f); end
  if isempty(findstr(fname,'.eps')), fname=[fname '.eps']; end

  cdir = pwd;

% ----------------------------------------------------------------------
% cd /home/weichsel/public_html/ML

% if first few letters of fname match the name of subdir, go there
  if 0
     d = dir('./');
     ii = find(cat(2,d.isdir)); d = d(ii);
 
     d2 = [];
     for i=1:length(d)
         if d(i).name(1)~='.'
             if strncmpi(fname, d(i).name, length(d(i).name))
                d2 = d(i).name;
             end
         end
     end
     if ~isempty(d2), eval(['cd ' d2]), end
  end
  disp([ 10 pwd 10])

% keyboard

  if isempty(fileparts(fname)), fname=['./' fname]; end
  if exist(fname)
     inp = input (['File ' fname ' exists. Overwrite? <[1]|0> ']);
     if isempty(inp); inp=1; end
     if inp~=1; disp ' '; return; end
  end

% disp([10 'Print current figure (' num2str(gcf) ') to file ' fname ' ...' ])
  fprintf(1,['\n   Print current figure (%g,%s) to ' fname ...
  '\n   pwd : %s\n\n'], gcf, get(gcf,'Renderer'), pwd);

% print -depsc2 -tiff myfile.eps - see info under __info.m
  CMD = sprintf( ...
  'print(gcf, ''-depsc2'', ''-tiff'', ''%s'')', fname);

% ptype = get (gcf, 'PaperType'); % save old setting (reset again below)
% set (gcf, 'PaperUnits', 'normalized');
% set (gcf, 'PaperType', 'A5');
% set (gcf, 'PaperSize', ... ); % NB! PaperSize is read only!
  set (gcf, 'PaperPositionMode','auto')% WYSIWYG (get onscreen figure size)

  h=findall(groot,'type','axes','tag','frame');;
  if ~isempty(h)
     set(h,'visible','off') % takes off frame line, yet keeps headers
  end

  disp([CMD 10])
  eval(CMD);

% if ~nodisp
  eval(['! gv ' fname ' &']); % end

  eval(['cd ' cdir]);
% set(gcf,'PaperType', ptype ) % set back to original size

  if ~isempty(h), set(h,'visible','on'); end

end

