function [s,q]=num2rat(xx,fmt,sep)
% Function:[s,i] = num2rat(x [,fmt, sep])
%
%    print number in rational format using wbrat.mex
%    if rational approximation fails, float number is
%    printed using fmt (%.5g)
%
% Wb,Dec29,09 ; Wb,Jan12,11

% for simple version based on matlab's rat
% see Archive/num2rat_110112.m

  if nargin<2 || isempty(fmt), fmt='%.4g'; end
  if nargin<3 || isempty(sep)
     if isvector(xx), sep=', '; else sep=' '; end
  end

  xx(find(abs(xx)<1E-12))=0;
  [z,q]=wbrat(xx,'-q');

  ss=cell(size(xx)); n=numel(xx);
  for i=1:n, x=xx(i);

     if q.ee(i)==0
        if q.Q(i)==1
             ss{i}=sprintf(fmt,q.P(i));
        else ss{i}=sprintf([fmt '/%g'],q.P(i),q.Q(i)); end
     else
        [z(i),q1]=wbrat(x^2,'-q');
        if q1.ee==0
           if x<0, s='-'; else s=' '; end
           if q1.Q==1
                ss{i}=[s sprintf('sqrt(%g)',q1.P)];
           else ss{i}=[s sprintf('sqrt(%g/%g)',q1.P,q1.Q)];
           end
        else ss{i}=sprintf(fmt,x); end
     end
   end

   for i=1:n, ss{i}=[sep ss{i}]; end

   if isvector(ss)
      s=cat(2,ss{:}); s=s(numel(sep)+1:end);
   else
      s=zeros(size(ss));
      for i=1:n, s(i)=numel(ss{i}); end;  n=max(s(:));
    % for i=1:n, ss{i}=[repmat(' ',1,n-s(i)), ss{i}]; end
      for i=1:n, ss{i}(end+1:n)=' '; end
      for i=1:size(ss,1), ss{i,1}=cat(2,ss{i,:}); end
      s=cat(1,ss{:,1});
   end

end

