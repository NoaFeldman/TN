function [s,i]=system_getline(cmd)
% function s=system_getline(cmd)
% Wb,Aug21,12

  [i,s]=system(cmd);

% NB! system (as well as evalc('! cmd') capture intermediate
% input from the MatLab prompt @#(*&! (eg. by typing command
% while previous command is still running => skip!
% affects: repHome, pwd, cd, ...
% Wb,Aug21,12

  l=max(find(s(1:end-1)==10)); % allow one trailing newline
  if ~isempty(l), s=s(l+1:end); end

end

