function sh=mp(varargin)
% Function sh=mp(M [,OPTS]))
%
%    surface plot of matrix (if matrix is 3D, plot layer by layer)
%    where M is a 2D or 3D matrix; returns surface handle sh.
%
% Explicitely specify x-, y-, z-data using 'xd' and 'yd' below
% or by the following alternative usgage:
%
%    sh=mp([xd [,yd [,zd]],] M [,OPTS]))
%
% Options
%
%     'nsecs', nsecs    number of seconds to wait from one plot to the next
%                       <0: pause | =0: drawnow (default) | >0: wait nsecs
%     'xla', xlabel     specify label for x-axis
%     'yla', ylabel     specify label for y-axis
%     'zla', zlabel     specify label for z-axis
%     'iz', zidx        index of z direction to use (default: all)
%     'zlim', [z1 z2]   ZLim for surface plot (default: 'auto')
%     'clim', [c1 c2]   CLim for surface plot (default: 'auto')
%     'view', [al be]   view for surface plot (default: 2D view)
%     'diag', val       set diagonal to val (eg NAN or 0, ...)
%     '3D'              3D view
%     'smooth'          3D view
%     'cm'              set colormap cmap0w
%     'xy'              use axis xy instead of axis ij
%     'trans'           transpose (ij) data
%     'lh'              call light2() (-> cmap irrelevant!)
%     'addrc'           adds 1 column and row (for surface() to show all data)
%     'trim', <nr>      trim <nr> of cols/rows around initial matrix
%     'tstr', ...       title string
%
%     '-log[xyz]'       (semi)log[xyz], i.e. '-logx','-logy', and '-logz'
%     '-gca'            plot into current axes set (same as mp0)
%
% Wb,Jun24,08 ; Wb,Nov26,02 ; (C) AW 2002

% DEPRECATED: ------------------------------------------------------- %
%     'xd', xdata       specify data for x-axis (see leading args)
%     'yd', ydata       specify data for y-axis
%     'zd', ydata       specify data for z-axis
% ------------------------------------------------------------------- %

% check whether xdat, ydat or zdat is specified as leading arguments
  for n=1:nargin, if ~isnumeric(varargin{n}), n=n-1; break, end; end

  if nargin<1 || isempty(n) || n<1 || n>4
     eval(['help ' mfilename]);
     dispstack(dbstack) % NB! may happen through recursive call of mp(!)
     if nargin || nargout, error('Wb:ERR','invalid usage (%g,%g)',n,nargin)
     else wblog('ERR','invalid usage (%g/%g)',n,nargin); 
     end, return
  end

  M=varargin{n};
  if n>1, ydat=varargin{1}; else ydat=1:size(M,1); end
  if n>2, xdat=varargin{2}; else xdat=1:size(M,2); end
  if n>3, ydat=varargin{3}; else zdat=1:size(M,3); end
  varargin=varargin(n+1:end);

  getopt('init',varargin);
     gcaflag=getopt('-gca');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if gcaflag, hold off; else
  if ~strcmp(get(gca,'Tag'), 'MP_AXES') || narg>0
  %| length(get(gca,'UserData'))<2
     na = narg;
     first_call = 1;

     clf; gca;

     [n1 n2 n3] = size(M);
     
     if n3==1 && ~isreal(M)
         M(:,:,3) = M(:,:,1).*conj(M(:,:,1));
         M(:,:,2) = imag(M(:,:,1));
         M(:,:,1) = real(M(:,:,1));
         n3 = 3;
         
     end
     if n3>1

        slh = makeslider([1, n3, 1/(n3-1), 20/(n3-1)]);

        set(gca, 'Tag', 'MP_AXES');
        set(gca, 'UserData', {double(M), slh, varargin{:}});

        set(slh, 'Callback', [get(slh,'Callback') 'mp; ']);

       % colormap button

         cm = uicontrol('Style', 'togglebutton', 'String', 'CM', 'Tag', 'mpCM', ...
             'Units', 'norm', 'Position', [0.02 0.02 0.08 0.08]);
		
         set(cm, 'UserData', get(gcf, 'Colormap'));
         set(cm, 'Callback', [ ...
         'if get(gco,''Value''); cmap0w; else ' ...
         'set(gcf, ''Colormap'',get(gco,''UserData'')); ' ...
         'set(gca, ''CLimMode'',''auto''); end; ']);
		
		% lighting button
		
         lb = uicontrol('Style', 'togglebutton', 'String','Light', 'Tag','mpLT',...
             'Units', 'norm', 'Position', [0.02 0.10 0.08 0.08]);
		
         set(lb, 'Callback', [ ...
         'if get(gco,''Value''); ' ...
         '   set(gco, ''UserData'', get(gcf,''Colormap'')); light2; ' ...
         'else ' ...
         '   delete(findall(gca,''Type'', ''light'')); ' ...
         '   set(gcf, ''Colormap'',get(gco,''UserData'')); ' ...
         '   set(gca, ''CLimMode'',''auto''); ' ...
         'end; ']);
		
		% 2D button
		
         v2 = uicontrol('Style', 'togglebutton', 'String', '2D', ...
             'Units', 'norm', 'Position', [0.02 0.30 0.08 0.08]);
		
         set(v2, 'Callback', [ ...
         'if get(gco,''Value''); ' ...
         '   [az el] = view; set(gco, ''UserData'', [az el]); view(2); ' ...
         'else; view(get(gco,''UserData'')); end; ']);
		
		% 3D button
		
         v3 = uicontrol('Style', 'togglebutton', 'String', '3D', ...
             'Units', 'norm', 'Position', [0.02 0.22 0.08 0.08]);
		
         set(v3, 'Callback', [ ...
         'if get(gco,''Value''); ' ...
         '   [az el] = view; set(gco, ''UserData'', [az el]); view(3); ' ...
         'else; view(get(gco,''UserData'')); end; ']);
     end
  else
     ud = get(gca, 'UserData');

     M   = ud{1};
     slh = ud{2};
     varargin=ud(3:end);

     na = length(ud)-1;
     first_call = 0;
  end

  end % gcaflag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ndims(M)<2, help mp, dispstack(dbstack)
      return
  end

  [n1 n2 n3] = size(M);

  getopt ('init', varargin);
     nsecs  = getopt ('nsecs', -1   );

     zidx   = getopt ('zidx', 1:n3); % old
     zidx   = getopt ('iz',   zidx);

     xdat   = getopt ('xd', xdat);
     ydat   = getopt ('yd', ydat);
     zdat   = getopt ('zd', zdat);

     xla    = getopt ('xla','X');
     yla    = getopt ('yla','Y');
     zla    = getopt ('zla','Z');

     logx   = getopt ('-logx');
     logy   = getopt ('-logy');
     logz   = getopt ('-logz');

     zlimx  = getopt ('zlim',  []); % (NB! zlim is a MatLab function!)
     climx  = getopt ('clim',  []); %
     viewx  = getopt ('view',  []); % (NB! view is a MatLab function!)
     sdval  = getopt ('diag',  []);
     trim   = getopt ('trim',   0);
     tstr   = getopt ('tstr','z =');

     trans  = getopt ('trans'     );
     cmap   = getopt ({'cmap','cm'});
     f3D    = getopt ('3D'        );
     fsmooth= getopt ('smooth'    );
     xyflag = getopt ('xy'        );

     lght   = getopt ('lh'        );

     addrcf = getopt ('addrc');
     addrcf=1; % addrc used most of the times anyway (Wb,Aug20,08)
     if getopt ('~addrc'), addrcf=0; end

  if (getopt('check_error')), return, end

  if strcmpi (zlimx, 'auto'), zlimx=[]; end
  if strcmpi (climx, 'auto'), climx=[]; end
  if strcmpi (viewx, 'auto'), viewx=[]; end

  if ~isempty(zlimx), if length(zlimx)~=2
     wblog('WRN - invalid ZLim parameter (dim=%d)', length(zlimx));
     zlimx = [];
  end, end

  if ~cmap, 
     h = findall(gcf, 'Type', 'uicontrol', 'Tag', 'mpCM');
     if ~isempty(h), if get(h(1),'Value'), cmap=1; end; end
  end

  if ~lght, 
     h = findall(gcf, 'Type', 'uicontrol', 'Tag', 'mpLT');
     if ~isempty(h), if get(h(1),'Value'), lght=1; end; end
  end

% same for clim (instead of zlim) ...
  if ~isempty(climx),  if length(climx) ~=2
      wblog('WRN - invalid CLim parameter (dim=%d)', length(climx));
      climx = [];
  end, end

  if cmap && ~isempty(climx)
      wblog('WRN - either clim or cmap (ignore later)');
      climx = [];
  end

  M=double(M); if issparse(M), M=full(M); end

  if max(size(M))==1, addrcf=1; end

  if sum(intersect(zidx,1:n3)~=zidx)
     wblog('\nERR - Index out of range: %d (%d, %d)\n', zidx(1), 1, n3)
     return
  end

  ilast = get(get(gca, 'zlabel'), 'UserData');

  if ~gcaflag && first_call && exist('slh')
     i = zidx(1);
     set(slh, 'Value', i);

   % store ilast in gca->zlabel->userdata
     set(get(gca, 'zlabel'), 'UserData', i);

     set(gcf,'CurrentObject',slh);
     eval(get(slh,'Callback')); %% calls mp again where it proceeds!

     return

  elseif length(zidx)>1
   % just accept zidx only!
   % --------------------------------------------------------
     is = get(slh, 'Value');

     if is==ilast, i = ilast;
     elseif is>ilast
          [ii,i] = findmel(zidx(find(zidx>ilast)), is, -1); % closest
     else [ii,i] = findmel(zidx(find(zidx<ilast)), is, -1); % closest
     end
     if isempty(i), i = ilast; else i = i(1); end

     if i~=is
        set(slh, 'Value', i);

        set(gcf, 'CurrentObject', slh);
        eval(get(slh,'Callback')); %% calls mp again where it proceeds!

        return
     end
  else
     if exist('slh')
          i = get(slh,'Value');
     else i=1; end
  end

  set(get(gca, 'zlabel'), 'UserData', i);
% wblog('i=%d, is=%d, ilast=%d',i,is,ilast);

  M2 = M(:,:,i);

       if ~isempty(sdval), M2 = M2-diag(diag(M2)+sdval);              end
       if trim>0,          M2 = M2(trim+1:end-trim, trim+1:end-trim); end
       if trans,           M2 = M2';                                  end

  if findstr(lower(get(gca, 'NextPlot')), 'replace'), cla; end
  if addrcf
      if length(xdat)>1
           xdat = [xdat(1)-diff(xdat(1:2)), mrow(xdat)];
      else xdat = [xdat-1, xdat]; end

      if length(ydat)>1
           ydat = [ydat(1)-diff(ydat(1:2)), mrow(ydat)];
      else ydat = [ydat-1, ydat]; end

      sh  = surface( xdat, ydat, addrc(M2));
  else sh = surface( xdat, ydat,       M2 ); end
  shading flat

  if f3D, view(3); end
  if fsmooth, shading interp; end
  
  if ~isempty(climx), set(gca, 'CLim', climx); end
  if ~isempty(zlimx), set(gca, 'ZLim', zlimx); end
  if ~isempty(viewx), view(viewx); end

% axis equal
  axis ij, xytight,  box on
  if lght, light2;
  elseif cmap, cmap0w; end

  if xyflag, axis xy, end
    
  mmax = max(max(M2));
  mmin = min(min(M2));

  lims = [mmin, mmax, mmax-mmin];
  fact = max(floor(log10(setdiff(abs(lims),0))));

  if abs(fact)>3
     lims=lims/(10^fact);
     tstr=sprintf('%s [%.2g %.2g; %.2g]*10^{%d}',tstr,lims(1:3),fact);
  else
     tstr=sprintf('%s [%.2g %.2g; %.2g]',tstr,lims(1:3));
  end

  if length(zdat)>1 || ~strcmp(upper(zla), 'Z')
  tstr=sprintf('%s=%g (%s)', zla, zdat(i), tstr); end

  label3(xla,yla,'',tstr);
% wblog('k=%d --> i=%d', k, i);

  if logx, set(gca,'XScale','log'); end
  if logx, set(gca,'XScale','log'); end
  if logz, set(gca,'ZScale','log'); end

% view(3); % view([-26 8]);
% set(gca,'ZLim', [-.1 12]); xytight

if nargout<1, clear sh; end

% figure(gcf); NB! redraws entire figure

return

