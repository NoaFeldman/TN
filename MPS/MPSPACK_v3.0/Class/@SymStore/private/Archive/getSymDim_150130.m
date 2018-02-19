function dd=getSymDim(varargin)
% function dd=getSymDim(varargin)
%    
%    usage: getSymDim <copy/paste string of CData/CGRef info>
% => use regexp to extract symmetry type and labels of input
%
% Wb,Jan16,15

  i=-1; for i=1:nargin, if ~ischar(varargin{i}), i=-1; break; end, end
  if i<0
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  istr=[varargin{:}]; t=istr; vflag=0;


  t=regexprep(t,'.*(S[pU])\(*(\d+)\)\s*','$1$2\n');

  t=regexprep(t,'\[([+-]+);[^\]]*\]\s*','$1\n');

  t=regexprep(t,'(\([\d,]+\))[,\s]*','$1\n');

  t=regexprep(t,'w=.*','');
  t=strread(t,'%s','whitespace','\n')';

  D=[getenv('RC_STORE') '/' t{1} '/RStore'];
  if ~exist(D,'dir')
     error('Wb:ERR','\n   ERR RCStore/%s does not yet exist',t{1}); end

  qs=strread(t{end},'%s','whitespace','(,; )\t\n');

  if numel(t)<1, istr, t
     error('Wb:ERR','\n   ERR invalid usage'); end

  if numel(t)==1
     ff=dir(D);
       ff=ff(find([ff.isdir]==0));
       [~,is]=sort([ff.bytes]); ff=ff(is);
     nf=numel(ff); qs=cell(1,nf); vflag=2;

     for i=1:nf
        qs{i}=regexprep(ff(i).name,'\((\d+)\).*','$1');
     end
  end

  if ~nargout && vflag<=0, vflag=1; end
  if vflag, fprintf(1,'\n'); end

  nq=numel(qs); dd=zeros(1,nq);
  for i=1:nq
     mat=sprintf('%s/(%s).rep',D,qs{i});
     if ~exist(mat,'file'), error('Wb:ERR',...
        '\n   ERR multiplet %s (%s) not found',t{1},qs{i});
     end
     q=load(mat,'-mat'); q=q.RSet;

     dd(i)=size(q.Z,1); t=regexprep(q.type,'[\(\)]','');
     if vflag
        if all(q.J>=0 & q.J<=9 & q.J==round(q.J))
             qstr=sprintf('%g',q.J);
        else qstr=vec2str(q.J,'-f'); end

        istr=sprintf('%s  (%s)  d=%d',t,qstr,dd(i));

        if vflag>1, istr={istr,datestr(ff(i).datenum)};
               fprintf(1,'  %3d.  %-24s   %s\n',i,istr{:});
        elseif fprintf(1,'  %3d.  %s\n',i,istr);
        end
     end
  end

  if ~nargout, clear dd; fprintf(1,'\n'); end

end

