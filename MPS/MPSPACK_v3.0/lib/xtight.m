function xl=xtight(varargin)
% xtight - set x axis tight to data
% usage: xtight([fac, opts])
%
%    fac = extra zoom (default fac=1, fac>1 introduces white margin)
%
% Options
%
%   'view'   zoom to xlim within local ylim view
%   'x1',..  fixed xlim(1)
%   'x2',..  fixed xlim(2)
%
% Wb,2002

  if nargin && isnumeric(varargin{1})
  fac=varargin{1}; varargin=varargin(2:end); else fac=[]; end

  getopt ('init', varargin);
     x1 = getopt('x1',[]);
     x2 = getopt('x2',[]);
     vw = getopt('view' );
     dflag= getopt('-d' );
  varargin=getopt('get_remaining');

  ah=gca;

  mh=[ findall(ah,'tag','xmark'); findall(ah,'tag','ymark') ];
  if ~isempty(mh)
   % set temporarily to invisible
     setprops(mh,'-s','Visible','off');
  end

  if vw
     xl=getxlim('-view');
     if diff(xl)<=0
      % if xl(1)==xl(2)
      %    if xl(1), xl=xl(1)*[.5 2]; else xl=[-1 1]; end; return
      % else error('Wb:ERR','getxlim() returned [%s]',vec2str(xl)); end
        error('Wb:ERR','getxlim() returned [%s]',vec2str(xl));
     end
     xlim(xl)
  elseif dflag
     xl=getxlim('-data');
     set(ah,'XLim',xl);
  else
     xopts={'YLim', get(ah,'YLim'), 'YLimMode', get(ah,'YLimMode')};
     axis tight; set(ah,xopts{:}); xl=xlim;
  end

  if ~isempty(fac)

     if isequal(get(ah,'XScale'),'linear')
        dx=(fac-1)/2*diff(xl); set(ah,'XLim', [xl(1)-dx, xl(2)+dx ]);
     else
        if xl(1)==0, xl=getxlim('-data'); end
        if all(xl~=0)
           fac=exp((fac-1)/2*diff(log(abs(xl))));
           if fac~=0 && ~any(isinf(xl)) && ~(isnan(fac) || isinf(fac))
           set(ah,'XLim', [xl(1)/fac, xl(2)*fac]); end
        end
     end

  elseif isequal(get(ah,'YScale'),'log') && xl(1)==0
     set(ah,'XLim',getxlim('-data','-pos'));
  end

  if ~isempty(x1) || ~isempty(x2)
     xl=xlim;
     if ~isempty(x1), xl(1)=x1; end
     if ~isempty(x2), xl(2)=x2; end
     xlim(xl);
  end

  if ~isempty(mh), setprops(mh,'-reset'); end
  if nargout, xl=xlim; else clear xl; end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

