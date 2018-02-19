function s=hostname
% Usage: n=hostname()

  s=evalc('!hostname');
  s=strtrim(s); % skip newline

end

