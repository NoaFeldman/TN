function s=itags2str(t)
% function s=itags2str(t)
% Wb,Feb24,16

% outsourced from info.m

   if isempty(t), s='';
   elseif iscell(t), t=sprintf(', %s',t{:});
      s=['{ ' t(3:end) ' }'];
   elseif ischar(t),
      s=['{ ' regexprep(t,'[,;| ]+',', ') ' }'];
   else 
      t, error('Wb:ERR','\n   ERR invalid itags');
   end

end


