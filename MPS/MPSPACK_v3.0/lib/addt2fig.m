function  h=addt2fig(varargin)
% Function h=addt2fig([tstamp, ] sflag)
% add time stamp at the right top corner of gcf
% default time is now(), but can be altered by numeric value tstamp.
%
%     default: Dec11,07 11:43:58
%     'M'      Muenchen, Mon dd, 2007
%     'H'      hostname, Mon dd, 2007
%     'Wb'     AW, Mon dd, 2007
%     'MWb'    AW, Muenchen, Mon dd, 2007
%
% Wb,Dec01,00

  if nargin && isnumeric(varargin{1})
  t=varargin(1); varargin=varargin(2:end); else t={}; end
  narg=length(varargin);

  if narg>0, sflag=varargin{1};
     tstr=time(t{:},'-l');
     switch(sflag)
        case 'M'
           tstr=['Muenchen, ', tstr(1:12)];
        case 'H'
           tstr=[getenv('HOST'), ', ', tstr(1:12)]; tstr(1)=upper(tstr(1));
        case 'MWb'
           tstr=['AW, Muenchen, ', tstr(1:12)];
        case 'Wb'
           tstr=['AW, ', tstr(1:12)];
     otherwise
        wblog('ERR','\NInvalid flag %s\N', sflag);
        error(' ');
     end
  else
     tstr=time(t{:});
  end

  h = header('hright', tstr); % time string at the right corner of the figure
  set (h(1), 'FontSize', 8);

  if nargout==0, clear h, end

end

