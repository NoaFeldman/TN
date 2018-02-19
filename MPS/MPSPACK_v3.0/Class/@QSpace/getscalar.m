function aa=getscalar(A)
% Function a=getscalar(A)
% Wb,Aug28,08

  aa=zeros(size(A));
  for k=1:numel(A), Ak=A(k);
     if numel(Ak.data)==1 && numel(Ak.data{1})==1
        if ~isempty(Ak.Q) && size(Ak.Q{1},1)>1 % safeguard
        wblog('WRN','severe QSpace inconsistency !??'); end

        a=Ak.data{1};

      % NB! rely on CData/CStore normalization convention
      % that |cdata|^2=1 => need to include cgr.cgw only
        if ~isempty(Ak.info) && isfield(Ak.info,'cgr') && ~isempty(Ak.info.cgr)
           cgr=Ak.info.cgr; 
           for i=1:numel(cgr), cgr(i).cgw=norm(mpfr2dec(cgr(i).cgw)); end
           a=a*prod(cat(1,cgr.cgw));
        end

     elseif (isempty(Ak.Q) || isempty(Ak.Q{1})) && ...
            (isempty(Ak.data) || isempty(Ak.data{1}))
        a=0;

     else error('Wb:ERR','not a scalar !??'); end

     q=[abs(a), normQS(A(k))]; % safeguard // Wb,May08,13
     if any(q~=0)
        e=abs(diff(q)/max(q)); if e>1E-8
           error('Wb:ERR','\n   ERR severe inconsistency in scalar value!');
        end
     end

     aa(k)=a;
  end

end

