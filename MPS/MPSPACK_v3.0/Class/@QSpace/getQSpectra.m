function [q,EE,qq,dz,DZ]=getQSpectra(H,ws)
% function [q,EE,qq,dz,DZ]=getQSpectra(H,ws)
%
%   get symmetry-resolved energy spectra.
%
% Input
%
%   H  diagonal scalar operator (with diagonals stored only)
%   ws specifies which symmetry to pick
%
% Output
%
%   q   unique set of labels for symmetry ws
%   EE  cell array of energies, one cell for each symmetry in q
%   qq  full set of symmetry labels (expanded to same size as EE)
%   dz  full set of multiplet degeneracies
%   DZ  full set of multiplet degeneracies (expanded to same size as EE)
%
% Wb,Apr12,12

   Q=H.Q; if numel(Q)~=2 || ~isequal(Q{:})
      error('Wb:ERR','\n   ERR invalid (scalar) operator'); end
   if isdiag(H)~=2, error('Wb:ERR',['\n   ' ...
     'ERR got non-diagonal input operator !??']); end

   if isreal(ws), ws=getsym(H,'-I',ws); end

   DZ=getzdim(H);
   DZ=expand2cell(DZ(:,:,1),H.data);
   QQ=expand2cell(Q{1},H.data);

   H=struct(H);

 % EE: eigenspectrum resolved w.r.t. specified symmetry sector only
   Q=round(1E6*H.Q{1})*1E-6; n=size(Q,1);
   [q,I,D]=uniquerows(Q(:,ws.j));

   m=numel(I); EE=cell(m,1); dz=cell(m,1); qq=cell(m,1);

   for j=1:m
      [EE{j},is]=sort(cat(2,H.data{I{j}}));
      dz{j}=cat(1,DZ{I{j}}); dz{j}=dz{j}(is,:);
      qq{j}=cat(1,QQ{I{j}}); qq{j}=qq{j}(is,:);
   end
 % EE=cat2(1,EE{:},{nan});

 % Wb,Sep07,15
   DZ=dz;
   for i=1:numel(dz)
      e=norm(diff(dz{i}));
      if e>1E-12, wblog('WRN',...
         'got varying dz within symmetry sector !? (e=%.3g)',e);
       % NB! e.g for the same abelian charge sector,
       % may have several different spin sectors! // Wb,Apr20,16
       % keyboard
      end
      dz{i}=mean(dz{i});
   end

   dz=[dz{:}]';

 % round dz values close to integer to integer
   i=find(abs(dz-round(dz))<1E-12); dz(i)=round(dz(i));

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function x=expand2cell(x,data)

 % x_=x;
   x=mat2cell(x,ones(size(x,1),1),size(x,2));
   for i=1:numel(data)
      x{i}=repmat(x{i},numel(data{i}),1);
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

