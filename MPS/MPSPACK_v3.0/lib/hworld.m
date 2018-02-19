function hworld(varargin)
% function hworld(varargin)
% Wb,Jan25,11

   fprintf(1,'\n  Hello world ...\n\n');

   for i=1:nargin
       fprintf(1,'  %4d %s\n',i,varargin{i});
   end

   fprintf(1,'ismcc: %g\n',ismcc);
   fprintf(1,'isdeployed: %g\n',isdeployed);
   fprintf(1,'usejava(awt): %g\n',usejava('awt'));
   fprintf(1,'usejava(jvm): %g\n',usejava('jvm'));
   fprintf(1,'isbatch: %g\n',isbatch);

end

