function s=add2struct(s,varargin)
% Function: s=add2struct(s,var1,var2,..)
%
%   Adds variables to given structure with the fieldname
%   being the variable name (therefore no expressions are
%   allowed as field values.
%
%   NB! input s may also be an emtpy object, or equivalently '-'.
%
% Variable name may be specified explicitely as
%
%   'var1'   variable with name `var1' must exist
%   'var1?'  variable with name `var1' (only if it exists)
%   'name1:var1[?]'  using `name1' as fieldname for variable `var1'
%
% Wb,Jul09,07

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if ~isstruct(s) && ~isempty(s)
     if isequal(s,'-'), s=struct; else
     error('Wb:ERR','Invalid first argument (not a structure)'); end
  end

  global flag__ val__

  for i=1:length(varargin)
      n=inputname(i+1);
      if ~isempty(n)
         s=setfield(s,n,varargin{i});
      elseif ischar(varargin{i}), n=varargin{i};
         if n(end)~='?', opt=0; else n=n(1:end-1); opt=1; end

         l=find(n==':' | n=='=');
         if ~isempty(l), l=l(1);
            nx=n(l+1:end); n=n(1:l-1);
            if n(end)=='?', n=n(1:end-1); opt=1; end
         else nx=n; end

         cmd=sprintf([...
           'global flag__ val__; flag__=0; ' ...
           'try, val__=%s; flag__=1; catch; end '], nx);
         evalin('caller',cmd);

         if flag__
            s=setfield(s,n,val__);
         elseif ~opt
            wblog('ERR','invalid variable or expression ''%s''',nx);
            l=lasterror; inl(1), disp(l.message), inl(1)
          % hint: eg. variables may not exist
         end
      else 
      error('Wb:ERR','failed to assign data (arg #%d)',i+1); end
  end

  clear global flag__ val__

  n=inputname(1);
  if nargout==0 && ~isempty(n)
     assignin('caller',n,s);
     clear s
  end

end

