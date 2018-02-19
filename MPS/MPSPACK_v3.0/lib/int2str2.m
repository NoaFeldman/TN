function cc=int2str2(dd)
% function cc=int2str2(dd)
%
%    introduce integer comma separator for thousand
%
% Wb,Aug09,10

  cc=cell(size(dd)); nel=numel(dd);
  for i=1:nel
      if abs(dd(i)-round(dd(i)))>1E-12
         wblog('WRN','data appears to contain non-integer (%s)',s);
         cc{i}=sprintf('%f',dd(i)); continue;
      end

      s=sprintf('%.0f',dd(i)); n=length(s);
      if n>3
         m=mod(n,3); if m, m=3-m; s=[repmat(' ',1,m),s]; end
         s=reshape(s,3,[]); s=[repmat(',',1,size(s,2)); s];
         s=s(m+2:end);
      end

      cc{i}=s;
  end
  if numel(cc)==1, cc=cc{1}; end

end

