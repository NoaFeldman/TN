function s=num2str2(val,varargin)
% function s=num2str2(val, lim1,scale1,tag1, lim2,tag2,scale2, ...)
%
%    returns sprintf('%.3g %s',val/lim_i,tag1) if val<lim_i
%
% alternative usage(2): s=num2str2('-n',val1,str1,val2,str2, ...)
%
%    eg. num2str2('-n',1,'error',0,'warning',...)
%    returns: 1 error and no warnings (i.e. this checks plural s)
%    Wb,Mar10,11
%
% alternative usage(3): s=num2str2(val,eps)
%
%    e.g. num2str2(0.45093,0.009)
%    returns: 0.450(9)
%    to be read as: the last digit '0' comes with an uncertainty +/- '9'
%
% Wb,Apr17,09

% tags: bytes2str, byte2str, sizestr

  getopt('init',varargin);
     vflag=getopt('-v');
     fmt  =getopt('fmt',[]);
     bflag=getopt('-bytes');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg==0 && numel(val)==1 && ~bflag
     if isempty(fmt), fmt='%.4g'; end
     if isreal(val)
        s=sprintf(fmt,val);
        if isempty(findstr(s,'.'))
           s=regexprep(s,'[eE][+0]*','E');
        else
           q=regexprep(s,'.*\.\d[eE][+0]*([1-9]+)','$1');
           if ~isequal(q,s), q=str2num(q)-1;
              s=sprintf([fmt 'E%g'],val*10^(-q),q);
           end
        end
     else
        r=real(val); i=imag(val);
        if r==0
             s=sprintf([fmt 'i'],i);
        else s=sprintf([fmt, regexprep(fmt,'%\+*','%+') 'i'], r,i);
        end
     end
     return
  end

% -------------------------------------------------------------------- %
% usage2
% -------------------------------------------------------------------- %

  if isequal(val,'-n'), nset=1;
  elseif isequal(val,'-n0'), nset=2; else nset=0; end

  if nset
     if mod(narg,2), error('Wb:ERR','invalid usage (%g)',narg); end
     s={};
     for i=1:2:narg
        q=varargin{i}; t=varargin{i+1};
        if numel(q)~=1 || ~isnumeric(q), error('Wb:ERR',...
           'invalid usage (alternating numbers expected)'); end
        if ~ischar(t), error('Wb:ERR',...
           'invalid usage (alternating string expected)'); end
        if q~=0 || nset>1
           if q==0 s{end+1}=sprintf('no %ss',t);
           elseif q==1, s{end+1}=sprintf('%g %s',q,t);
           else s{end+1}=sprintf('%g %ss',q,t); end
        end
     end
     if numel(s)>1
        s(2,:)={', '}; s{2,end-1}=' and '; s{2,end}='';
     end

     s=cat(2,s{:}); return
  end

% -------------------------------------------------------------------- %
% usage 3
% show number limited by given (absolute) accuracy
% -------------------------------------------------------------------- %

  if ~bflag && narg && isnumeric(varargin{1}) && isscalar(varargin{1})
     dval=varargin{1};

     if dval<0, error('Wb:ERR','invalid error margin (negative !??)'); end
     if dval==0, s=sprintf('%g(0)',val); return; end

     x=floor(log10(abs(dval/val)));
     if x>0 % uncerainty comparable or larger than value itself
     s=sprintf('%.3g +/- %.3g',val,dval); return; end

     if 0
      % show deviation of last significant digit
        if isempty(fmt)
           if abs(x)>4, fmt=sprintf('%%.%dg',abs(x)+2);
           else fmt='%g'; end
        end

        s=sprintf(fmt,val);
        ss=strread(s,'%s','whitespace','eE');

        if length(ss)>1
           i=find(s=='e'); if isempty(i), e='E'; else e='e'; end
           dval=dval/str2num(['1E' ss{2}]);
           ss={ss{1}, e, ss{2}};
        else
           ss={ss{1}, '', ''};
        end

        n=length(ss{1});
        i=findstr(ss{1},'.'); if isempty(i), i=n+1; end
        i=i-floor(log10(abs(dval)));

        if i>n
           ss{1}=[ss{1} '(0)'];
        else
           s=num2str(round(dval * 10^(-floor(log10(dval))) ));
           ss{1}=[ss{1}(1:i) '(' s ')'];
        end

      % fprintf(1,'%g @ %g\n',val,dval);
        s=cat(2,ss{:});
     elseif dval>abs(val)
         if dval>5*abs(val), s=sprintf('+/-%.2g',dval);
         else, s=sprintf('%.1g',val); s=['(' s(1) ')' s(2:end)];
         end
     else
      % show one digit beyond last significant digit
        if isempty(fmt)
           if abs(x)>4, fmt=sprintf('%%.%dg',abs(x)+2);
           else fmt='%.16g'; end
        end
        ss={ sprintf(fmt,val)
             sprintf(fmt,val+dval)
             sprintf(fmt,val-dval) };
        s=sum(diff(double(strvcat(ss)),1).^2,1);
        i=regexp(ss{1},'[eE]'); n=numel(s);
        if ~isempty(i)
           if numel(i)~=1 || any(s(i:end))
            % error('Wb:ERR','got no matching exponent!');
            % do not make programs stop because of this! // Wb,Aug27,15
              wblog('ERR','got no matching exponent !?');
              s=sprintf('[%g %g]',val,dval);
              return;
           end
        else i=n+1;
        end

        j=find(s,1); % first non-matching character
        if isempty(j) || j>=i
         % error('Wb:ERR','severe error within function');
           wblog('ERR','failed to represent accuracy of number');
           s=sprintf('<ERR %g!?>',val); return
        end

        if vflag
           b=repmat(' ',1,n); b(j)='|'; if i<=n, b(i)=':'; end
           printf('\n'); printf('  %s\n',ss{:},b);
           printf('\n');
        end

        s=ss{1};
        s=[s(1:j-1) '(' s(j) ')' s(i:end)];
     end
     return
  end

% -------------------------------------------------------------------- %
% default usage
% -------------------------------------------------------------------- %

  if bflag
     if ~narg
        varargin={ 
           5E3, 1,'bytes', ...
           1E6, 2^10,'kB', ...
           1E9, 2^20,'MB', ...
           1E12, 2^30,'GB', ...
           1E15, 2^40,'TB'
        };
        narg=numel(varargin);
     else error('Wb:ERR','invalid usage'); end
  end

  if narg<3 || mod(narg,3) || ~isscalar(val) || ~isnumeric(val)
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  varargin=reshape(varargin,3,[]); n=size(varargin,2);
  for i=1:n
     if ~isscalar(varargin{1,i}) || ~isnumeric(varargin{1,i}) || ...
        ~isscalar(varargin{2,i}) || ~isnumeric(varargin{2,i}) || ...
        ~ischar(varargin{3,i})
        error('Wb:ERR','invalid usage');
     end
  end

  [x,i]=sort(cat(2,varargin{1,:})); varargin=varargin(:,i);

  if isempty(fmt), fmt='%.4g %s';
  elseif isempty(regexp(fmt,'%.*s')), fmt=[ fmt ' %s']; end

  for i=1:n
     if val<varargin{1,i} || i==n
        s=sprintf(fmt,val/varargin{2,i},varargin{3,i});
        break;
     end
  end

end

