function s=num2tex(x,varargin)
% Function: num2tex(x [,'%..e'])
%
%    print number in %e format (%.4G)
%    replacing the E+NN part by 10^{NN}
%
% Wb,Mar07,13 : extended to vectors and matrices
% Wb,Oct10,05

  getopt('init',varargin);
     sep=getopt('sep',{', ', '\n'});
  fmt=getopt('get_last','%.4G');

  n=numel(x); s=cell(size(x));
  for i=1:n
      s{i}=num2tex_1(x(i),fmt);
  end

  if ~iscell(sep), sep={sep}; end
  if ~isempty(findstr(sep{1},'\')), sep{1}=sprintf(sep{1}); end

  if n==1, s=s{1};
  elseif n>1 && numel(find(size(x)>1))==1 % got vector
     s=join_to_string(s,sep{1});
  else % got matrix
     if numel(sep)==1, sep={sep{:},'\n'}; end
     if ~isempty(findstr(sep{2},'\')), sep{2}=sprintf(sep{2}); end
     for i=1:size(x)
        s{i}=join_to_string(s(i,:),sep{1});
     end
     s=join_to_string(s(:,1),sep{2});
  end

end

function s=join_to_string(s,sep)
  s=reshape(s,1,[]);
  s(2,:)={sep}; s{end}='';
  s=cat(2,s{:});
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function s=num2tex_1(x,fmt)

  if x==0, s=sprintf('%g',x); return; end

% for test purposes
  if ~ischar(x)
     i=find(fmt=='e'); % use E to enforce exponentials
     if ~isempty(i) && abs(log10(abs(x)))<3
        fmt(i)='g';
     end
     s=sprintf(fmt,x);
  else s=x; end

  s=regexprep(upper(s),'E[+-]00',''); % skip 10^0

  hasE=any(find(s=='E'));
  hasdot=any(find(s=='.'));

  if hasdot
   % skip trailing zeros right before and after an E, if there is an E
     s=regexprep(s,'[.0]*([E+-])0*','$1');
   % skip trailing zeros if there is no E (these are after the . !)
     if ~hasE, s=regexprep(s,'[.0]*$',''); end
  end

  if hasE
   % replace scientific represenation by 10^..
     s=regexprep(s,'E\+0*([0-9]*)','\\cdot10^{$1}');  % E+ (+ always present)
     s=regexprep(s,'E\-0*([0-9]*)','\\cdot10^{-$1}'); % - appears in exponent
     s=regexprep(s,'^1\.0*\\cdot','');  % 1.00E-12
     s=regexprep(s,'^1\\cdot','');      % 1E-12
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

