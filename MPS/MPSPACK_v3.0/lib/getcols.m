function ncols=getcols()
% function ncols=getcols()
% Wb,Nov11,11
   
% NB! after repeated calls, this may take 10 sec for each call !(*Y_!
% Wb,May02,13
  persistent nc
  if isempty(nc)
     nc=system('getcols.pl -q 2>/dev/null');
  end

  ncols=nc;

end

