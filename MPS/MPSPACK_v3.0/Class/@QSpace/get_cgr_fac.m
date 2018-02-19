function [dfac,s]=get_cgr_fac(A,i,sflag)
% function [dfac,s]=get_cgr_fac(A,i)
%
%    CGData is no longer stored with the QSpaces themselves,
%    only their weights cgw within outer multiplicity space.
%    Therefore CGData is globally normalized to 1. up to
%    outer multiplicity.
%
% Options
%
%    sflag=='s' => also return range as string info
%
% Wb,Oct17,14

% adapted from fromer get_cgr_range.m

  if nargin<2 || nargin>3 || numel(A)~=1 || numel(i)~=1
  error('Wb:ERR','\n   ERR invalid usage'); end

  if ~gotCGS(A), q=1; % [1;1];
  else
     if ~isfield(A.info,'cgr') || ~isfield(A.info.cgr,'size')
     error('Wb:ERR','\n   ERR invalid info.cgr data'); end

     cgr=A.info.cgr(i,:); q=ones(1,numel(cgr));
     for j=1:numel(cgr)
      % s=cgr(j).size;
      % d=cgr(j).qdir; r=numel(d);
      % d=(d=='+');
      % if sum(d)<=r/2, k=find(d); else k=find(~d); end
      % q(1,j)=prod(s(k));

        w=mpfr2dec(cgr(j).cgw);
        if isempty(w)
           if ~isempty(cgr(j).type) || ~isempty(cgr(j).qset) || ...
              ~isempty(cgr(j).qdir), error('Wb:ERR',...
              '\n   ERR unexpected CGRef data (assuming CGR_ABELIAN)');
           end
           w=1;
        end

      % check deviation from default normalization
        q(j)=norm(w); % Wb,Nov04,14
 
      % aw=abs(w); iw=find(aw==max(aw));
      % check deviation from default normalization sqrt(max(s))
      % e.g. see QSpace tag EXT_CGD_NORM // Wb,Oct28,14
      % q(j)=w(iw)/sqrt(max(s));
     end
  end

% dfac=1/sqrt(prod(q(1,:)));
% dfac is used by display.m => don't apply cgc factor on data sector
% as this may lead to confusion when actually looking up A.data{i};
% NB! also adapted NormCGC such that rank-2 CGC's are always stored
% as identity matrix (still), yet rank>2 CGC's have |cgc|^2 = 1
% (or delta_i,j in outer multiplicity i and j). % Wb,Oct18,14
  dfac=1;

  if nargin>=3 && ~isempty(sflag) || nargout>1
   % q=q(2,:);
     if norm(q-1)>1E-12, q=prod(q);
        if numel(q)>1
           if norm(diff(q))>1E-12
                s=sprintf(' * %.4g',q); s=s(4:end);
           else s=sprintf('%.4g (x%d)',q(1),numel(q)); end
        else
           s=sprintf('%.4g',q);
        end
     else s='';
     end
     if nargin<2, q=s; end
  end

end

