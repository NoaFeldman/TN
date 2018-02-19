function i=isempty(A)
% Wb,May12,06;  Wb,Aug23,08

  i=1; n=length(A);
  if ~builtin('isempty',A)
     for k=1:n, d=A(k).data; q=A(k).Q; id=1; iq=1; e=0;

        if ~isempty(d)
           if ~iscell(d), e=1; elseif ~isempty(d{1}), id=0; end
        end
        if ~isempty(q)
           if ~iscell(q), e=1; elseif ~isempty(q{1}), iq=0; end
        end

        if e || xor(id,iq) && ~isscalar(A)
           if n==1
                error('Wb:ERR','\nsevere QSpace inconsistency !??');
           else error('Wb:ERR','\nsevere QSpace inconsistency (%d/%g) !??',k,n);
        end, end

      % scalar may have iq=1, id=0
        if ~id, i=0; break; end
     end
  end

end

