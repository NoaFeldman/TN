function A=makeIrop(A)
% function A=makeIrop(A)
%
%    reattach 3rd irop dimension if not yet present
%    (this concerns only rank-2 scalar operators).
%
% Wb,Feb12,16

% see also: appendScalarSymmetry.m

  for k=1:numel(A)
     r=numel(A(k).Q); if ~r, continue; end
     if r==2
      % NB! do not just fix Q{} data since this misses
      % to properly also adjust cgr, itags, etc.!
        E3=getIdentityQS(A(k),getvac(A(k))); % s,0;s*
      % NB! got OUTgoing operator index, and make it last!
        A(k)=QSpace(contractQS(A(k),2,E3,'3*'));
     elseif r~=3
        error('Wb:ERR','\n   ERR invalid usage (got rank-%g op !?)',r);
     end
  end

end

