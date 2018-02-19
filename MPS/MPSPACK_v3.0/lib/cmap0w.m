function cmap0w(varargin)
% cmap0w - set colormap ((zdat=0) = white)
% NB! uniformly spaced similar region in colormap for neg AND pos values!
%
%     usage: cmap0w (options)
%
% Options:
%
%     'dim'    - dimension of resulting colormap
%     'cfac'   - factor on clim -> cfac<1 increases the contrast!
%     'fact'   - linear scaling for fact=1 (otherwise: (0:dc:1) .^ fact)
%     'axis'   - specify axis to look at for surface handles
%     'sh'     - surface handle explicitely specified
%     'all'    - set all handles to the same clim value (scale all the same)
%
% Wb, Oct12,02

    getopt ('init', varargin);
       ndim = getopt ('dim',  256);
       fact = getopt ('fact', 1. );
       sall = getopt ('all'      );
       cfac = getopt ('cfac', 1. );
       ah   = getopt ('axis', [] );
       sh   = getopt ('sh',   [] );
    varargin=getopt('get_remaining');

    o={}; if length(varargin)==1
       o={'-data', varargin{:}};
    elseif length(varargin)>1
       eval(['help ' mfilename]);
       error('Wb:ERR','invalid usage');
    end

    if ~isempty(sh),ah=get(sh,'parent'); end

    if cfac<=0
       eval(['help ' mfilename]); return
       error('Wb:ERR','invalid cfac = %g', cfac);
    end

    if isempty(ah)
       if sall
            ah=findall(gcf,'Type','axes')';
       else ah=gca; end
    end

    if isempty(sh), sh = [
       findall(ah,'Type','surface')
       findall(ah,'Type','patch')
    ]'; end

    cp = [0 0 1];  % positive (`b')
    c0 = [1 1 1];
    cn = [1 0 0];  % negative (`a')

    ndim = ndim - rem(ndim,2);  %% make ndim even!

    dc = 1/(ndim/2); %% ndim/2 * 2 for pos & neg region!
    xx = (0:dc:1) .^ fact;

    cm = [ ...
       ones(ndim/2+1,1) * c0 + flipud(xx'       ) * (cn-c0); ...
       ones(ndim/2  ,1) * c0 +       (xx(2:end)') * (cp-c0)];
    colormap (cm);

    if sall
        if isempty(sh)
           disp('ERR - no suface handles found.')
           eval(['help ' mfilename]); return
        end

        [a,b,cl]=getDataRange_(o{:},sh);

        if ~isempty(cl)
        set(h,'clim',cl*cfac); end
    else
        for h=ah
            sh = [ findall(h,'Type','surface'); findall(h,'Type','patch')]';

            [a,b,cl]=getDataRange_(o{:},sh);

            if ~isempty(cl)
            set(h,'clim',cl*cfac); end
        end
    end

  % keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% '-data' option useful for mp() for example when cmap should affect
% the whole 3D data block the same way

function [a,b,cl]=getDataRange_(varargin)

  if nargin>1 && isequal(varargin{1},'-data')
   % assume numeric input
     a=0; b=0;
     for i=2:length(varargin);
         dd = varargin{i}(:);
         a = min ([a, min(dd)]); % a<=0
         b = max ([b, max(dd)]); % b>=0
     end
  else
     a=0; b=0; sh=varargin{1};
     for h2=sh
         dd = get(h2,'CData');
         a = min ([a, min(min(dd))]); % therefore a<=0!
         b = max ([b, max(max(dd))]); % therefore b>=0!
     end
  end

  ab = max(abs(double([a b])));

  if ab~=0
       cl=[-ab, +ab];
  else cl=[]; end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

