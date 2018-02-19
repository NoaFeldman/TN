function xx=getylim(varargin)
% Function xx=getylim([OPTS])
% Options
%
%    '-view'   determine YLim fully shows all y-data in given view
%    '-data'   default: ylim is full y-data range
%
% Wb,Jan10,08 - see also ytight.m Wb,2002.

  if nargin>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(varargin), varargin{1}='-data'; end

  switch varargin{1}
    case '-view', xx=getYLimView;
    case '-data', xx=getYLimData;
    otherwise
       eval(['help ' mfilename]);
       error('Wb:ERR','invalid usage')
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% determine ylim only within xlim view
% skip line handles with viewer than 3 data points
% since quite likely these are just markers

function yy=getYLimView()

  lh=findall(gca,'Type','Line','visible','on');
  m=length(lh); if ~m, yy=[]; return; end

  yy=nan(m,3); xl=xlim;

  for i=1:m
      xd=get(lh(i),'XData'); ii=find(xd>=xl(1) & xd<=xl(2));
      if isempty(ii), continue; end

      yd=get(lh(i),'YData'); yd=yd(ii);
      yy(i,:)=[min(yd), max(yd), length(xd)];
  end

% NB! if only straigth lines with only few data points exist
% this may lead to somewhat unexpected behavior!

  i=find(yy(:,3)>2); if isempty(i), i=1:size(yy,1); end
  yy=[ min(yy(i,1)), max(yy(i,2)) ];

  if yy(1)==1E99
     wblog('WRN','failed to determine YLim for data within view');
     yy=ylim;
  end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% NB! on log scale ylim(1) might be =0 even though data won't be plotted !!

function yy=getYLimData()

  lh=findall(gca,'Type','Line','visible','on');
  m=length(lh); yy=zeros(m,2);

  for i=1:m
      yd=get(lh(i),'YData'); yd=yd(find(yd>0));
      if ~isempty(yd)
           yy(i,:)=[min(yd), max(yd)];
      else yy(i,:)=[+inf -inf];
      end
  end

  yy=[ min(yy(:,1)), max(yy(:,2)) ];

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

