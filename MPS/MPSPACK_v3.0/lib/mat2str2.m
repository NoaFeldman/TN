function str = mat2str2(M, varargin)
% function s = mat2str2(M [,opts])
%
%    writes matrix M as string
%
% Options
%
%   'fmt',..   is the format string ('%10g')
%   'sep',..   is the separator string (' ')
%   'rowsep',..is the string to separate rows ('\n')
%   'istr',..  info/intro string
%   '-f'       no shortcuts (enforce full mode)
%   'notiny'   no tiny numbers on the numerical noise level
%   'phase'    abs|phase instead of real+imag
%   'nofac'    do not use overall factors pulled to the front (1 for scalar, 0 otherwise)
%
% See also existing MatLab routine mat2str().
% Wb,Jul11,03

  if nargin<1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  str=[];

  getopt('init',varargin);
     fmt    = getopt('fmt','%8g');
     sep    = getopt('sep',' ');
     rowsep = getopt('rowsep','\n');
     istr   = getopt('istr','');
     notiny = getopt('notiny');
     pflag  = getopt('phase');
     fflag  = getopt('-f');
     nofac  = getopt('nofac');
  getopt('check_error');

  s=size(M); r=numel(s); n1=s(1); n2=s(2); n=numel(M);
  if ~nofac & isscalar(M), if isreal(M), nofac=1; end; end

  if n==1
     str=sprintf(fmt,M); return
  end
  if ~fflag
     if all(diff(M(:))==0)
        str=sprintf([fmt ' (%s)'],M(1),vec2str(s,'sep','x')); return
     end
  end

  if r>2, error('Wb:ERR',...
    '\n   ERR invalid usage (got rank-%g object)',r); end
  if ~fflag && n1==n2
     d=diag(M);
     if norm(M-diag(d))==0
        if all(diff(d)==0)
             str=sprintf([fmt ' (eye; %gx%g)'],M(1),s);
        else str=sprintf('diag([%s])',vec2str(d,'fmt',fmt)); end
        return
     end
  end

% eps = 5E3 * abs(2-sqrt(2)^2);
  eps = 1E6 * abs(2-sqrt(2)^2) * abs(max(M(:)));

% wblog('eps = %g', eps);
  
% get 10^x prefactor if numbers are too large/small
% fact = mean(abs(M(:)));
  if nofac
     fact=1;
  else
     fact = max (abs(M(:)));
     if fact~=0, fact = floor(log10(fact)); end
     if findstr(fmt, 'd'), fact=1; end
     if abs(fact) > 3
        fact = 10^fact;
        M = M / fact;
     else
        fact = 1;
     end
  end

% keyboard
   
  str=[];

  if isreal(M)
     for i=1:n1
       for j=1:n2
          if j==1
             if i==1
                  str = [ str sprintf([       fmt], M(i,j)) ];
             else str = [ str sprintf([rowsep fmt], M(i,j)) ]; end
          else    str = [ str sprintf([   sep fmt], M(i,j)) ]; end
       end
     end
  else                             % for example:  %+12.5g
      [fm1,r] =  strtok(fmt, '%.gefGEF');       %   +ِ12
      [fm2,r] =  strtok(r,   '%.gefGEF');       %       5
      [fmc,r] =  strtok(fmt, '%+-.0123456789'); %        g

      pdot = findstr(fmt,'.');
      if ~isempty(fm1) & ~isempty(pdot)
          if pdot<findstr(fmt,fm1)
             fm2 = fm1;
             fm1 = '';
          end
      end

    % wblog('fmt tokens: %s, %s, %s', fm1, fm2, fmc);

      if  ~isempty(fm1), fm1 = str2num(fm1); else fm1=8; end
    % else eval(['help ' mfilename]); return; end

      if ~isempty(fm2)
          fm2  = str2num(fm2);
          fmtr = sprintf('%%.%dg', fm2);
          fmtc = sprintf('%%+.%dg', fm2);
      else
          fm2 = [];
          fmtr = '%g';
          fmtc = '%+g';
      end
	  fmts = sprintf('%%%ds', fm1);

    % wblog('fmtr >%s<, fmts >%s<', fmtr, fmts);

      for i=1:n1
        for j=1:n2

          mij = M(i,j);
          if notiny
             if abs(real(mij))<eps, mij = imag(mij); end
             if abs(imag(mij))<eps, mij = real(mij); end
          end

          if real(mij)==0
             if imag(mij)==0
                vstr = '0.';
             else
                vstr = sprintf([fmtr 'i'], imag(mij));
             end
          else
             if imag(mij)==0
                vstr = sprintf(fmtr, real(mij));
             else
                if ~pflag
                vstr = sprintf([fmtr fmtc 'i'], real(mij), imag(mij));
                else
                vstr = sprintf([fmtr '|' fmtr], abs(mij), angle(mij)/pi);
                end
             end
          end

          if j==1
             if i==1
                  str = [ str sprintf([       fmts], vstr) ];
             else str = [ str sprintf([rowsep fmts], vstr) ]; end
          else    str = [ str sprintf([   sep fmts], vstr) ]; end
       end
     end
  end

% wblog('>%s<', rowsep);
% wblog('>%s<',    sep);

  if ~isempty(istr) && isempty(find(istr=='='))
  istr=[istr ' = ']; end

  if fact~=1 | ~isempty(istr)
     if fact~=1, vstr = sprintf('%1.0E * ', fact);
     else        vstr = ''; end

     if size(M,1)>1 & (~isempty(findstr(rowsep,'\n')) | ~isempty(findstr(rowsep,10)))
         if     ~isempty(find(istr=='[')), bs='\n]';
         elseif ~isempty(find(istr=='{')), bs='\n}'; else bs=''; end
         str = sprintf(['%s%s\n\n%s' bs], istr, vstr, str);
     else
         str = sprintf('%s%s[%s]',   istr, vstr, str);
     end
  end

return

