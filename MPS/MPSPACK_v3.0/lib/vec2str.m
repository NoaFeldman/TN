function str=vec2str(v, varargin)
% function str=vec2str(v [,opts])
%
%    write numerica vector as string
%
% Options
%
%    'fmt',...  number format (%g)
%    'sep',..   separation of values (' ')
%    '-f'       show full vector without any shortcuts such as x1:dx:x2
%
% Wb,Nov12,99, Jul08,03, Jan16,06

  str='';

  getopt ('init', varargin);
     fmt    = getopt ('fmt','%g');
     sep    = getopt ('sep', ' ');
     fflag  = getopt ('-f'); % formerly 'full' // Wb,Jul07,11
     nofac  = getopt ('nofac'   );
  if (getopt ('check_error')), return, end

  v=full(v);
  if ~isvector(v)
      if numel(v)==0, str='[]'; return; end
    % NB! mat2str is a standard MatLab routine!
      o={'fmt',fmt,'sep',sep}; if fflag, o{end+1}='-f'; end
      str=mat2str2(v,o{:});
      return
  end

  n=numel(v); if ~n, return; end
  if n==1, str=sprintf(fmt,v); return; end
  if n>1 && ~fflag && ~nofac && all(diff(v)==0)
      str=sprintf([ fmt ' (x%g)'],v(1),length(v));
      return
  end

% check if elements in v are equally spaced
  dxconst_flag = 0;

  eps = 1E-4;

  if fflag==0
  if n > 4
     diff_v = diff(v);
     avg_dv = mean(diff_v);
     std_dv = std (diff_v);

     if std_dv==0
        if avg_dv(1)==1
             str = sprintf([fmt ':' fmt],v(1), v(end));
        else str = sprintf([fmt ':' fmt ':' fmt],v(1), diff_v(1), v(end));
        end
        return
     end

     if std_dv/(abs(avg_dv)+eps) < eps  % & std_dv~=0
        dxconst_flag = 1;
     end
  elseif n==1
     str = sprintf(fmt,v(1));
     return
  end
  end

% get 10^x prefactor if numbers are too large/small
  fac = mean(abs(v));
  if fac~=0, fac = floor(log10(fac)); end
  if ~nofac && abs(fac) > 3
     fac = 10^fac;
     v = v / fac;
  else
     fac = 1;
  end

% keyboard
   
  if dxconst_flag
    
     delta = (max(v)-min(v)) / (n-1);
     if delta~=0
        dbl = 10^(round(log10(delta))+5);
        delta = (delta+dbl)-dbl;
     end

     str = [sprintf(fmt,min(v)) ':'];

     if delta == 1
        str = [ str sprintf(fmt,max(v)) ];
     else
        str = [ str sprintf(fmt,delta ) ':' sprintf(fmt,max(v)) ];
     end
  else
     str = [ sprintf(fmt,v(1)), sprintf([sep fmt],v(2:n)) ];
  end

  if fac~=1
     str = sprintf('10^{%.0g} * [%s]', log10(fac), str);
  end

end

