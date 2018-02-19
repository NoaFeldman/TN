function [ff,I]=getnrgfiles(nrg)
% function [ff,I]=getnrgfiles([nrg])
%
%    get filenames of NRG file structure (default nrg: './NRG/NRG')
%    if requested, I returns data of info file.
%
% Wb,Aug22,12

  if nargin>1 || nargin==1 && ~ischar(nrg)
      error('Wb:ERR','invalid NRGdata tag'); end
  if ~nargin, nrg='./NRG/NRG'; end

  if ~isempty(findstr(getenv('USER'),'Weichsel')) && nargout<2
     ff=dir2([nrg '_\d+.mat']);
     return
  end

  ff=dir([ nrg '_*.mat']); n=numel(ff); N=[]; I=[];

  i=find(nrg=='/');
  if ~isempty(i), i=i(end);
       ftag=nrg(i+1:end); p=nrg(1:i);
  else ftag=nrg; p=[];
  end

  if ~n, error('Wb:ERR','invalid NRGdata tag'); end

  fmt='%s_%02d.mat';
  if isempty(findstr(ff(1).name, [ftag '_00.mat']))
     if ~isempty(findstr(ff(1).name, [ftag '_000.mat']))
          fmt='%s_%03d.mat';
     else error('Wb:ERR','invalid NRGdata tag');
     end
  end

% check file order
  for i=1:n
     s=sprintf(fmt,ftag,i-1);
     if ~isempty(findstr(ff(i).name,s))
       % NB! matlab removes path in field ff(i).name
         if ~isempty(p), ff(i).name = [p, ff(i).name]; end
     else break; end
  end
  N=i-1; if ~N
     disp(strvcat(ff(i:end-1).name)); error('Wb:ERR',['\n   ' ...
    'ERR invalid NRGdata (%s: missing files)'],nrg);
  end

% make sure, info file exists
  for i=i:n
     if ~isempty(findstr(ff(i).name,[ftag '_info.mat']))
        I=ff(i).name; if ~isempty(p), I = [p, I]; end
        break
     end
  end
  if isempty(I)
     disp(strvcat(ff(i:end-1).name)); error('Wb:ERR',['\n   ' ...
    'ERR invalid NRGdata (%s: missing info file)'],nrg);
  end

  if nargout>1, I=load(I); end

  ff=ff(1:N);

end

