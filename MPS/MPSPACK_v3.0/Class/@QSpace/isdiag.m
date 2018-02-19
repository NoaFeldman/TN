function [isd,estr]=isdiag(A,dflag)
% function [isd,estr]=isdiag(A)
%
%    check whether QSpace is in diagonal rank-2 format.
%    return value:
%       1   if block-diagonal in symmetry space
%       2   if, in addition, all data{} sets are vectors, e.g.
%           compressed format as returned bei QSpace::EigenSymmetric()
%       0   otherwise
%
%    If all data blocks are of dimension 1x1, by default isd=1
%    is returned, unless dflag=='-d' is specfied, which returns isd=2.
%
% Wb,Apr24,09

  n=length(A.data); isd=-1; estr='';

  if length(A.Q)~=2, estr=sprintf(...
     '%s requires rank-2 object (%d).',mfilename,length(A.Q));
  elseif ~isequal(A.Q{:}), estr=sprintf(... 
     '%s requires block-diagonal operator.',mfilename); end
  if nargout<2 && ~isempty(estr), error('Wb:ERR',estr); end

  for i=1:n, ai=A.data{i}; s=size(ai);
     if length(s)~=2 || s(1)~=s(2) && all(s~=1)
        estr=sprintf('expecting symmetric operator!');
        break;
     end
     if isempty(ai), continue; end % Wb,Sep26,15

     switch isd
       case 2 % vectorized diagonal representation
         if all(s~=1), estr=sprintf(...
           'confusing block dimensions (compressed !??)'); break; end
       case 1
         if s(1)~=s(2)
            if any(s==1), estr=sprintf(...
              'confusing block dimensions (compressed !??)'); break
            else, estr=sprintf(...
              'got rectangular data block (%s)',vec2str(s)); break
            end
         else e=norm(ai-diag(diag(ai)));
            if e>1E-12
             % estr=sprintf('got non-diagonal data (%.3g)',e);
               isd=0; break
            end
         end
       case 0
       otherwise
         if any(s>1)
            if s(1)==s(2), isd=1;
            elseif any(s==1), isd=2;
            else, isd=0; estr=sprintf(...
              'got rectangular data block (%s)',vec2str(s)); break
            end
       % otherwise can't say
         end
     end

     if isd==1
      % e=norm(ai-ai'); // Wb,Apr24,14
        e=norm(ai-diag(diag(ai)));
        if e>1E-12, isd=0; end
     end
  end

  if ~isempty(estr)
     if nargout<2, error('Wb:ERR',estr); end
     isd=0;
  elseif isd<0 % i.e. got all 1x1 data blocks
     if nargin>1
        if ~isequal(dflag,'-d')
           error('Wb:ERR','\n   ERR invalid dflag');
        end
        isd=2;
     else
        isd=1;
     end
  end

end

