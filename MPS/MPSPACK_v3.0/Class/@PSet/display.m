function display(P)

   if P.n==0
      fprintf(1,'\n   Empty parameter set.\n\n');
      return
   end

   fprintf(1,'\n  PARAMETER LENGTH  DATA\n\n');
   for i=1:P.r
      if length(P.data{i})<8,
      s=['[ ' mat2str2(P.data{i},'fmt','%8.5g') ' ]'];
      else s=sprintf('(%s array)',vec2str(size(P.data{i}),'sep','x')); end
      fprintf(1,'  %8s %5g    %s\n', P.name{i},P.s(i), s);
   end
   fprintf(1,'\n');

end

