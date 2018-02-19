function [wd,TT]=getWeightsNRG(varargin)
% function [wd,TT]=getWeightsNRG([nrg,]TT)
%
%    Get FDM weight distribution for given NRG run
%    (default: $LMA/NRG/NRG.mat).
%
% Wb,Jun20,13

  getopt('init',varargin);
     vflag=getopt('-v');
     TT=getopt('T',-1);
  varargin=getopt('get_remaining'); narg=length(varargin);

  if numel(varargin) && ischar(varargin{1})
       nrg=varargin{1}; varargin=varargin(2:end);
  else nrg=[getenv('LMA') '/NRG/NRG']; end

  if numel(varargin)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  f=sprintf('%s_info.mat',nrg);
  if ~exist(f,'file'), error('Wb:ERR',...
    '\n   ERR invalid NRG data (inaccessible info file) !??');
  end

  Inrg=load(f);

  if vflag, o={}; else o={'-q'}; end
  [ED,IN]=getEDdata(nrg,'-d',o{:}); N=numel(ED); nT=numel(TT);

  dloc=IN.dloc;
  if dloc<2 || dloc~=round(dloc)
     error('Wb:ERR','\n   ERR invalid d_loc = %g !??',dloc);
  end

  for k1=1:N
     if ~isempty(ED(k1).E), break; end
  end

  E0=min(cat(1,ED.E)); if E0~=0
     wblog('WRN','got small energy E0=%.3g !??',E0);
     for k=k1:N, ED(k).E = ED(k).E - E0; end
  end

% the following is fast (even for many T's)
% as compared to getEDdata() above! // Wb,Jun20,13

  for it=1:nT, T=TT(it);
     if T<0, if T==-1, T=-10; end
        T=Inrg.Lambda^(-(N+T)/2);
        TT(it)=T;
     end
     rw=zeros(N,1);

     if T==0, rw(end)=1;
     else
        for k=k1:N, q=ED(k); if ~isempty(q.E)
           rw(k) = q.deg' * exp(-q.E/T + log(dloc)*(N-k));
        end, end
        rw=rw/sum(rw);
     end

     wd(:,it)=rw;
  end

% keyboard

end

