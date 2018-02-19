function [xr,I]=getrhoQS(H,varargin)
% Function getrhoQS(H [,beta,opts])
%
%    get density matrix with respect to H for given beta
%    relative to local NRG energy scale.
%
% Options
%
%    '-sqrt'  also return sqrt(R) as I.r in I-structure => e.g. trace(rr')=1
%
% Wb,Jun02,08

  getopt('init',varargin);
     rflag=getopt('-sqrt');
     qflag=getopt('-q');
     vflag=getopt('-v');
     reps =getopt('-nnz',1E-8);
  beta=getopt('get_last',[]);

  H_=H; eflag=0;

% check whether H needs to be diagonlized first
  for i=1:numel(H.data), q=H.data{i}; s=size(q);
     if numel(find(s>1))<=1, continue; end % got vector
     if any(diff(s)), s
        error('Wb:ERR','\n   ERR got rectangular matrix !?'); end
     if norm(q-diag(diag(q)))<1E-14, H.data{i}=diag(q)';
     else 
        if ~qflag, wblog(1,'WRN','diagonalizing H'); end
        [ee,Ie]=eigQS(H_); H=Ie.EK; eflag=1;
        break;
     end
  end

  H=QSpace(H); zdim=prod(getzdim(H),2);

  ee=H.data; ee=cat(2,ee{:});
  Emin=min(ee); Emax=max(ee); dE=Emax-Emin;

  if isempty(beta)
   % Wilson site k=0 may have dE0~1E-2 with Emin~=0
     if dE>1E-6, beta=1E5/dE; else beta=100; end
  end
  if beta<0, error('Wb:ERR','invalid beta (%g)',beta); end
  
  R=H; nrm=0;
  for i=1:length(H.data)
     R.data{i}=exp(-beta*(H.data{i}-Emin));
     nrm=nrm+zdim(i)*sum(R.data{i});
  end

  nrm=1/nrm; if rflag, r=R; end
  nnz=0;

  for i=1:length(R.data)
     q=nrm*R.data{i}; q(find(abs(q)<1E-16))=0;
     R.data{i}=diag(q);  if rflag
     r.data{i}=diag(sqrt(q)); end
     nnz=nnz+numel(find(q>reps));
  end

  if eflag
   % return to old basis
     R=QSpace(contractQS(contractQS(Ie.AK,2,R,1),2,Ie.AK,2));
  end

% xr=R;
  xr=skipzeros(R);

  I=add2struct('-',nnz);
  if rflag, I.r=skipzeros(r); end
  if qflag, I.msg=''; return; end

% check uniqueness of ground state
  d=diag(xr,'-d'); n=length(find(d>.5));
  if beta>10
    if n>1
       I.msg=sprintf('WRN: groundstate is degenerate (%g) !??',n);
       wblog('WRN',I.msg);
    else
       I.msg=sprintf('groundstate is unique (%g @ %g)',n,1-max(d));
    end
  else I.msg=''; end

end

