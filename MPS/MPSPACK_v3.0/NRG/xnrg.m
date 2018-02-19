
  U=0.12; epsd=-U/3; Gamma=0.01; B=0; Lambda=2;
  N=49; Nkeep=256;

  global param
  param=struct( ...
     'U', U, 'epsd', epsd, 'Gamma', Gamma, ...
     'B', B, 'N', N, 'Lambda', Lambda ...
  );

  if exist('BX','var'), param.BX=BX; end

  vflag=(N<=5); % store full unitary matrix ("MPS state")
  if exist('VFLAG','var'), vflag=VFLAG; end

% overall factor for consistency with standard NRG
% NB! this makes couplings converge to 1. within NRG iterations
  Hfac = 2/(1+Lambda);
   
% impurity Hamiltonian
% [HQ,H0] = getNRG_H0('fac', Hfac);

  H0 = diag([ 0, epsd+B/2, epsd-B/2, 2*epsd + U ]' + U/2);
  if isfield(param,'BX')
     H0(2,3)=param.BX/2;
     H0(3,2)=param.BX/2;
  end

  H0=H0*Hfac;

% couplings
% ff = getNRGcoupling(Gamma, Lambda, N,'-AL') * Hfac;
  ii=0:N-3; Delta=1/Lambda;
  xi = (1-Delta.^(ii+1))./sqrt((1-Delta.^(2*ii+1)).*(1-Delta.^(2*ii+3)));
  ff = [ sqrt(2*Gamma/pi), (1+Delta)/2*(Delta.^(ii/2)) .* xi ];

  ff = (ff*Hfac) .* (Lambda.^((1:length(ff))/2)); % rescale energies

  param  % display structure

% one particle state ops
  c2=[0 1; 0 0];          % the corresponding destruction operator
  e2=eye(2);

% z2=diag([1,-1]);        % |empty>=(1,0), |occ>=(0,1) comes second (!)
  z2=c2*c2' - c2'*c2;

% order: 0, up, down, up+down (i.e. up is fast index)
% tensor product (up,down) with up the fast index -> comes 2nd in kron()

% z4=diag([1 -1 -1  1]);
  z4=kron(z2,z2); e4=eye(4);

% fu = [ 0  1  0  0
%        0  0  0  0
%        0  0  0  1
%        0  0  0  0 ];
  fu = kron(e2,c2);

% fd = [ 0  0  1  0 
%        0  0  0 -1       % NB! up particle is added last, i.e.
%        0  0  0  0       % |up,down> = c^+_{up} c^+_{down} |0>
%        0  0  0  0 ];
  fd = kron(c2,z2);

  [vv,ee]=eig(H0); ee=diag(ee); E0=min(ee); ee=ee-E0;

  Fu = vv' * (ff(1)*fu) * vv;
  Fd = vv' * (ff(1)*fd) * vv;

  HN=cell(1,N); HN{1}=H0;
  EE=nan(256,N); EE(1:length(ee),1)=sort(ee);
  E0=[E0, zeros(1,N-1)];
 
% local operator set to remember
  if 1 || ~uflag, CN={ Fu Fd }; end

  CN=CN(ones(N,1),:);

  if vflag, VV=vv; end;

% Wb,Sep06,06: can use "tilde-operators" z4*f[ud] instead of just
% f[ud] on first site only (!!) cf. handnotes on calculating
% local operators opdates for correlation functions later on
  rL=sqrt(Lambda);

  fprintf(1,'\n');

  for k=2:N
      d=length(ee); Z=zeros(d); H=diag(rL*ee);

      H = kron(e4,H) ...
        + kron(z4*fu,Fu') + kron((z4*fu)',Fu) ...
        + kron(z4*fd,Fd') + kron((z4*fd)',Fd);

      % Hx= [ H   Fu' Fd' Z
      %       Fu  H   Z  +Fd'
      %       Fd  Z   H  -Fu'
      %       Z  +Fd -Fu  H 
      %     ];
      % if ~isequal(Hx,H), wblog('TST','debug ...'); keyboard, end

      HN{k}=H;
      fprintf(1,'   NRG %2d/%d: %d/%d ...\r',k,N,d,size(H,1));
            
      [vv, ee] = eig(H);
      [ee,idx] = sort(diag(ee)); E0(k)=ee(1); ee=ee-ee(1);
      
      m=min(length(ee),size(EE,1));
      EE(1:m,k)=ee(1:m);

      if length(ee)>Nkeep
       % keep degenerate blocks together
         for l=Nkeep+1:length(ee), if ee(l)>ee(Nkeep)+1E-8, break; end, end
         ee=ee(1:l-1);
         idx=idx(1:l-1);
      end

      vv=vv(:,idx);

      for s=1:size(CN,2)
      CN{k,s}=vv' * kron(eye(4),CN{k-1,s}) * vv; end

      if k<N
         E=ff(k)*eye(d);

         Fu = vv' * kron(fu,E) * vv;
         Fd = vv' * kron(fd,E) * vv;
      end

      if vflag
         VV=kron(e4,VV)*vv;
      end
  end

  fprintf(1,'\n\n');

  if ~exist('plotflag','var') || plotflag
     figure
     h1=plot(EE(:,1:2:end)','b'); hold on
     h2=plot(EE(:,2:2:end)','r');

     legend([h1(1),h2(1)],'odd','even','Location','NorthEast');

     xlabel('Wilson shell k');
     ylabel('rescaled eigenenergies');
     title('energy flow diagram');
  end

