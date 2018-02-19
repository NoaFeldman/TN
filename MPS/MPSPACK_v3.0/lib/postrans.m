function [pos,ta]=postrans(pos)
% Function [p,ta]=postrans(pos)
%
%   Translates position such as 'NorthEast', ... (similar to legend)
%
% Output
%
%   pos   position in normalized units
%   ta    text alignment (e.g. as used in postext)
%
% Wb,Aug15,08

  if nargin~=1, eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if iscell(pos)
     if numel(pos)~=2 || ~ischar(pos{1})
        error('Wb:ERR','invalid usage'); end
     z=pos{2}; pos=pos{1};
  else z=[]; end

  if ischar(pos)
     X=[0.05, 0.5, 0.95]; % default position offsets from axis set
     Y=[0.05, 0.5, 0.95];

     switch lower(pos)
       case {'north',    'n' }, i=[2,3];
       case {'northeast','ne'}, i=[3,3];
       case {'east',      'e'}, i=[3,2];
       case {'southeast','se'}, i=[3,1];
       case {'south',    's' }, i=[2,1];
       case {'southwest','sw'}, i=[1,1];
       case {'west',      'w'}, i=[1,2];
       case {'northwest','nw'}, i=[1,3];
       case {'center',   'cc'}, i=[2,2];

       otherwise error('Wb:ERR','invalid position tag `%s''',pos);
     end

     pos=[X(i(1)),Y(i(2))];

     if ~isempty(z)
        n=numel(z); if n>2, error('Wb:ERR','invalid usage (n=%g)',n); end
        if numel(z)==2, pos=pos+z;
        else
           if i(1)==1 || i(1)==3
                pos(2)=pos(2)+z;
           else pos(1)=pos(1)+z; end
        end
     end
  else
     if numel(pos)~=2, error('Wb:ERR','invalid position specification'); end
  end

  if nargout<2, return; end

  x=pos(1); y=pos(2);

  if x<0.25, ha='left';   elseif x>0.75, ha='right'; else ha='center'; end
  if y<0.25, va='bottom'; elseif y>0.75, va='top';   else va='middle'; end

  ta={'HorizontalAl',ha,'VerticalAl',va};

end

