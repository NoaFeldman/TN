function ah=ax2nfig(varargin)
% function ah=ax2nfig([ah0,][opts])
%
%    takes axis handle given (ah0; or current axis if not specified) and
%    copies the whole graph into a new figure;
%
% Options
%
%   -l   also copy legend axis set
%   -h   also copy header axis set,
%        while reducing the size of the actual axis set to be copied
%
% Wb,Apr03,01 ; Wb,Aug14,11

  if ~nargin || ~isnumber(varargin{1}) || ~isaxis(varargin{1}), ah0=gca;
  else
     ah0=varargin{1}; varargin=varargin(2:end);
  end

  getopt('init',varargin);
     hflag=getopt('-h');
     lflag=getopt('-l');
     zoom =getopt('zoom',0);
  getopt('check_error');

  set(gca,'selected','off')
  saveDisp(ah0);

  pos=get(groot,'DefaultAxesPosition');
  if hflag && ~zoom, zoom=0.9; end % pos([2 4])=[0.13 0.78];
  if zoom
     pos=[ pos(1:2)+((1-zoom)/2)*pos(3:4), zoom*pos(3:4) ];
  end

  f=figure;

  if hflag
     h=findall(get(ah0,'parent'),'tag','frame');
     if ~isempty(h)
        h2=copyobj(h,f);
     end
  end

  ah=copyobj(ah0,f);
  set(ah,'Units','Normalized','Position',pos);

  restoreDisp(ah0);
  restoreDisp(ah);

% copy also legend if available
  if lflag
    lh=findleg(ah0);
    if ~isempty(lh)
       lp=get(lh,'Position');
       p0=get(ah0,'Position'); l2=ax_pos2pos(lp,p0,pos);
       repax(axes('Position',l2), lh);
    end
  end

  if nargout==0, clear ah; end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function saveDisp(ah)
% save displayname to userdata (will not be copied otherwise!)
  for h=findall(ah,'type','line')'
   % NB! arrowsafe is a line object but has DisplayName disabled!
     if isfield(get(h),'DisplayName')
        setuser(h,'disp',get(h,'DisplayName'));
     end
  end
end

function restoreDisp(ah)
% restore displayname from userdata
  for h=findall(ah,'type','line')'
     if isfield(get(h),'DisplayName')
        set(h,'DisplayName',getuser(h,'disp','-rm'))
     end
  end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

