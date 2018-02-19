function f=fieldname(P,i)
% function f=fieldname(P,i)
% Wb,Nov22,11

  n=numel(P.name);
  if i>0 && i<=n, f=P.name{i};
  else
     wblog('ERR','field-id out of bounds (%g/%g)',i,n);
     f='';
  end

end

