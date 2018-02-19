function C=plus(A,B)
% overloading + operator
% Wb,Nov17,11

  if isa(A,'SymOp') && isa(B,'SymOp'), C=A;
     for i=1:numel(A)
        as=regexprep(A(i).istr,'^[^=]+=[ ]+','');
        bs=regexprep(B(i).istr,'^[^=]+=[ ]+','');

        na=nnz(A(i).op); nb=nnz(B(i).op);
        if na && nb, C(i).istr=[ as ' + ' bs ];
        elseif na, C(i).istr=as;
        elseif nb, C(i).istr=bs; end

        C(i).op=C(i).op+B(i).op;
        C(i).hc=C(i).hc+B(i).hc;
     end
  elseif isnumeric(B) && numel(B)==1 && isa(A,'SymOp')
   % treat scalar as identity operator

     C=A; q=C.op; C.istr=sprintf('%s %+g',C.istr,B);
     C.op=q+diag(repmat(B,1,size(q,1)));

  else error('Wb:ERR','\n   ERR invalid usage'); end

end

