function C=minus(A,B)
% overloading - operator
% Wb,Nov17,11

  if isa(A,'SymOp') && isa(B,'SymOp')
     as=regexprep(A.istr,'^[^=]+=[ ]+','');
     bs=regexprep(B.istr,'^[^=]+=[ ]+','');

     C=A; C.istr=[ as ' - ' bs ];
     C.op=C.op-B.op;
     C.hc=C.hc-B.hc;
  elseif isnumeric(B) && numel(B)==1 && isa(A,'SymOp')
     C=A; q=C.op; C.istr=['rmtrace( ' C.istr ' )'];
     C.op=q-diag(repmat(B,1,size(q,1)));
  elseif isnumeric(B) && isequal(size(A.op),size(B))
     as=regexprep(A.istr,'^[^=]+=[ ]+','');
     C=A; C.istr=[ as ' - <M>' ];
     C.op=A.op-B;
  else error('Wb:ERR','\n   ERR invalid usage'); end

end

