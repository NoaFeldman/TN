function cr=checkSpSz(Sp,Sz,type)
% function cr=checkSpSz(Sp,Sz,type)
% Wb,Jun28,10

  m=numel(Sz); n=numel(Sp); cr=zeros(m,n);

  for i=1:m
  for j=1:n
    % NB! Frobenius norm
    % in case of sparse object, where norm(S) results in an error!)
      q=comm(Sz(i),Sp(j)); if norm(q)<1E-12, continue; end
      cr(i,j)=olap(q,Sp(j)) / olap(Sp(j),Sp(j));

      e=norm(q-cr(i,j)*Sp(j)); if e>1E-12
      error('Wb:ERR','invalid Sp/Sz pair (%d/%d)',i,j); end
  end
  end

  if nargin<3, return; end

  if ~isempty(regexp(type,'^Z\d+$')) type='ZN'; end
  % see also is_abelian_symmetry // Wb,Aug29,16

  switch type

    case {'A','P','ZN'}
       if ~isempty(Sp) || numel(Sz)~=1 || ~isempty(cr)
       wblog('ERR','invalid Abliaan Sz/Sp setting'); end

    case 'SU2'
       if numel(Sz)~=1 || numel(Sp)~=1 || ...
        ~isnumber(cr) || abs(cr-1)>1E-12
       wblog('ERR','invalid SU(2) Sz/Sp pair (%.3g)',abs(cr-1)); end

    case {'SU3','SU4','SU5','SU6','SU7','SU8'}
     % added SU4-6: Wb,Apr04,14

       r=str2num(type(3:end))-1; % Wb,Apr13,12
     % NB! 3rd Sp can be extracted from existing two Sp
     % since: CP_13 = comm(CP_12,CP_23) already contained
       if numel(Sz)~=r || numel(Sp)~=r
       error('Wb:ERR','invalid %s operators',type); end

     % from initial test run
       switch r
          case 2
           % e=norm(cr.^2 - 0.25*[4 1; 0 3],'fro'); % with sqrt(3) factor in Z2
             e=norm(cr - [1 -0.5; 0 1.5],'fro');
          case {3,4,5,6,7,8,9} % Wb,Apr04,14
          % r=2: cr = [1 -0.5; 0 1.5];
          % r=3: cr = [1 -0.5 0; 0 1.5 -1; 0 0 2];
             e=(0:.5:r); q=diag(e(3:2+r)) + diag(-e(2:r),1);
             e=norm(cr - q,'fro');
          otherwise e=-1;
       end
       if abs(e)>1E-12, error('Wb:ERR','invalid %s operators',type); end
       
     % e=norm(comm(Sp(1),Sp(2))-Sp(3),'fro'); if e>1E-12, error(...
     % 'Wb:ERR','invalid SU3 operators: failed to generate 3rd Sp'); end
       for i=1:r
       for j=2:r
          Sij=comm(Sp(i),Sp(j)); e=trace(Sij);
          if abs(e)>1E-12, error('Wb:ERR',...
            'got trace-full %s generator! (%g,%g; %.3g)',type,i,j,e);
          end
       end
       end

    case {'Sp4','Sp6','Sp8','Sp10'} % Wb,Nov29,11 % Sp4: Wb,Apr13,12
       r=str2num(type(3:end))/2;

       if numel(Sz)~=r || numel(Sp)~=r
       error('Wb:ERR','invalid %s operators',type); end

     % from initial test run (cr(1:2,1:2) is identical with SU3 cr)
       if r==2
          e=norm(cr - [2 -2; 0 2],'fro');
       elseif r==3
        % e=norm(cr - [1 -0.5 0; 0 1.5 -2; 0 0 1],'fro');
          e=norm(cr - [2 -1 0; 0 3 -4; 0 0 2],'fro');
       elseif r==4
          e=norm(cr - [2 -1 0 0; 0 3 -2 0; 0 0 4 -6; 0 0 0 2],'fro');
       else e=-1; end
       if abs(e)>1E-12, error('Wb:ERR','invalid %s operators',type); end
       
       S3=comm(Sp(1),Sp(2)); e=trace(S3); if abs(e)>1E-12
       error('Wb:ERR','got trace-full 3rd %s generator! (%.3g)',type,e); end

    otherwise, error('Wb:ERR','invalid symmetry (%s)',type);
  end

end

