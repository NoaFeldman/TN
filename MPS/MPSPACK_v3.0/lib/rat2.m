function [str,e]=rat2(x,varargin)
% function str=rat2(x [,fmt])
% find fractional representation of x, if any, based on MatLab's rat.m
%
%    The format specifier fmt (default: %g) is used only
%    if no sensible fractional representation can be found.
%
% See also rats()
% Wb,Mar31,05

  getopt('init',varargin);
     udat =getopt('-u',{});
     sflag=getopt('-s');
     rflag=getopt('-r');
  fmt=getopt('get_last','%g');

  if numel(x)==1 && nargin==(1+sflag+rflag)  % Wb,Jun05,16
     if rflag
          str=sprintf('sqrt(%s)',regexprep(rats(x^2),' ',''));
     else str=regexprep(rats(x),' ','');
     end

     if sflag && x>0, str=['+' str]; elseif rflag && x<0 str=['-' str]; end
     if nargout>1, eval(['q=' str ';']); e=abs(q-x)/abs(x); end
     return
  end

  if nargin<1, eval(['help ' mfilename]); end
  if nargin<2, fmt='%g'; end

  if ~isempty(udat)
     if ~iscell(udat) || numel(udat)~=2
     error('Wb:ERR','invalid specification of units ''-u'',{uval,ustr}'); end
     if udat{1}>0, x=x/udat{1};
     else
      % wblog('WRN','got %g as unit - ignore',udat{1}); 
        udat={};
     end
  end

  if rflag
       [e,d]=rat(x^2);
  else [e,d]=rat(x);
  end

  if abs(e)<32 & abs(d)<32
     if e~=0
        if isempty(udat)
           if d~=1, str=sprintf('%d/%d',e,d);
           else str=sprintf('%d',e); end
        else 
           if e==1
              if d~=1, str=sprintf('%s/%d',udat{2},d);
              else str=sprintf('%s',udat{2}); end
           elseif e==-1
              if d~=1, str=sprintf('-%s/%d',udat{2},d);
              else str=sprintf('-%s',udat{2}); end
           else
              if d~=1, str=sprintf('(%d/%d)%s',e,d,udat{2});
              else str=sprintf('%d%s',e,udat{2}); end
           end
        end
     else str='0'; end
  else 
     if isempty(udat), str=sprintf(fmt,x);
     else str=[ sprintf(fmt,x) udat{2}]; end
  end

  if rflag, str=sprintf('sqrt(%s)',str); end
  if sflag && x<0, str=['-' str]; end
  if nargout>1, eval(['q=' str ';']); e=abs(q-x)/abs(x); end

end

