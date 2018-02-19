function ah=nrg_header(varargin)

  global paras param savemat

  evalin('caller','setglobal savemat');
% set global on workspace --> needed for addfinfo.m

  getopt('init',varargin);
     dy=getopt('dy',0);
     XP=getopt('-x','');
  getopt('check_error');

  try   checkParam(0);
  catch l % l=lasterror;
     wblog('ERR','checkParam: %s',l.message);
  end

  structexp(param)

  if 0 && isfield(param,'Lambda') && ...
     sum(isfield(param,{'J','JH','Gamma'}))==1 && ~isfield(param,'U3')
   % xor in J, JH, Gamma, ...

   % -------------------------------------------------------------- %
   % NB! passivate (missing extra/new fields in param) => take next %
   % -------------------------------------------------------------- %

     global gTK; nTK='T_K';
     if ~isempty(gTK)
        if iscell(gTK)
         % e.g. gTK = { 'T_K^{\chi}', TK0 }
           TK=gTK{2}; nTK=strrep(gTK{1},'\','\\');
        else TK=gTK; end
     else TK=TKondo; end

     sTK=sprintf([nTK '=%s'], num2tex(TK,'%.3E'));

     if ~exist('D','var')
        if exist('hbdim','var'), D=hbdim; end
     end

     if ~exist('L','var')
        if exist('Psi','var'), L=(length(Psi)-2)/2; end
        if exist('N','var'), L=N; end
     end

     if exist('Psi','var'), ss='DMRG'; else ss='NRG'; end

   % sprintf('\\epsilon_d/U=%s', rat2(epsd/U))
     if ~exist('D','var') && exist('Nkeep','var'), D=Nkeep; end

     lam={sprintf('\\Lambda=%g',Lambda)};
     if isfield(param,'z') && param.z, lam{end+1}=sprintf('z=%.3g',z); end
     if isfield(param,'ALambda') && ~isempty(param.ALambda)
        if ischar(param.ALambda), lam{end+1}='A_\Lambda';
        else lam{end+1}=sprintf('A_\Lambda=%g',ALambda); end
     end
     if numel(lam)>1, lam=[ lam{1},' (',strhcat(lam(2:end),'-s',', '),')'];
     else lam=lam{1}; end

     ss = {
        sprintf('%s  N=%g, %s, D=%d',ss,L,lam,D);
        ''
        sTK
      % U/\\pi\\Gamma=%.5g = U/pi/Gamma
     };

     if exist('U','var') && exist('Gamma','var')
        s={sprintf('\nSIAM: U=%g', U)};
        if exist('epsd','var')
        s{end+1}=sprintf(', \\epsilon_d=%s', vec2str(epsd, 'sep',', ')); end
        s{end+1}=sprintf(', \\Gamma=%s', mat2str2(Gamma,'fmt','%.5g','rowsep','; '));
        if exist('J','var'), s{end+1}=sprintf(', J=%g', J); end
        ss{2}=cat(2,s{:});
     elseif exist('J','var') || exist('fJ','var')
        if ~exist('J','var'), J=fJ;
        elseif exist('fJ','var') && J~=fJ
        wblog('ERR','J inconsistency (%g,%g)',J,fJ); end

        ss{2}=sprintf('; Kondo: J=%g',J);
     end

     if exist('Bfield','var') && ~exist('B','var'), B=Bfield; end
     if exist('B','var')
        if TK~=0 && ~isnan(TK)
             ss{2}=[ss{2}, sprintf(', B=%.3g TK', B/TK)];
        else ss{2}=[ss{2}, sprintf(', B=%.3g', B)]; end
     end
     ss=sprintf('%s %s (%s)', ss{:});

  elseif isfield(param,'Lambda') && any(isfield(param,{'J','JH','Gamma'}))

     xpat='nrgIO|wsys|ff|mat|sym|istr'; if ~isempty(XP), xpat=[xpat,'|',XP]; end
     ss=param2str(param,'-tex','-x',xpat);
   % ss=regexprep(ss,' ([a-zA-Z])([a-zA-Z]+)\s*=','$1_{$2}=');

  elseif isfield(param,'tLR') % Quantum measurement setup
     xpat='IO|mat|sym'; if ~isempty(XP), xpat=[xpat,'|',XP]; end
     ss=param2str(param,'-x',xpat);

     ss=regexprep(ss,'([^\\])Lambda','$1\\Lambda');
     ss=regexprep(ss,'([^\\])Gamma','$1\\Gamma');
     ss=regexprep(ss,'([^\\])delta','$1\\delta');
     ss=strrep(ss,'Vbias','V_{b}');
     ss=strrep(ss,'tLR','t_{LR}');
     ss=strrep(ss,'tdot','t_{dot}');
  elseif ~isempty(param)
     xpat='IO|mat|sym|istr|POP'; if ~isempty(XP), xpat=[xpat,'|',XP]; end
     ss=param2str(param,'-tex','-x',xpat);
  else ss=''; end

  ss=regexprep(ss,'(ALambda|A\\Lambda).*energies''','A_{\\Lambda}');
  ss=regexprep(ss,'Nkeep','N_K');
  ss=regexprep(ss,'Etrunc','E_{tr}');
  ss=regexprep(ss,'\<Nz\>','N_z');
  ss=regexprep(ss,'sym=','');

  ss=regexprep(ss,'([a-zA-Z])*=''[a-zA-Z]*(=|\\leq|\\geq)([^'']*)''','$1$2$3');
% if ~isempty(regexp(ss,'leq')), keyboard, end


  getbase Imain Inrg; % forward to / used by addfinfo() below
  if isempty(Imain) && isempty(Inrg)
   % try to get data from caller
     q0='---xvar-tmp---';

     setuser(0,'xvar',q0);
     evalin('caller','if exist(''Imain'',''var''), setuser(0,''xvar'',Imain); end');
     q=getuser(0,'xvar'); if ~isequal(q,q0)
        Imain=q;
        setuser(0,'xvar',q0);
     end

     evalin('caller','if exist(''Inrg'',''var''), setuser(0,''xvar'',Inrg); end');
     q=getuser(0,'xvar','-rm'); if ~isequal(q,q0)
        Inrg=q;
     end

   % whos Imain Inrg
  end

  if isfield(Imain,'finished')
     if ~isempty(findstr(Imain.finished,'not finished'))
        if exist('Idma','var') && isfield(Idma,'finished')
           Imain.finished=Idma.finished;
        elseif exist('Inrg','var') && isfield(Inrg,'finished')
           Imain.finished=Inrg.finished;
        end
     end
  end

  if isbatch, addt2fig Wb; else addt2fig; end
  fh=addfinfo;

  if isfield(param,'istr') && ~isempty(param.istr)
     s=dbstack; 
     if numel(s)>1, s=[ strrep(s(2).name,'_','\_') ': ']; else s=''; end
     [hh,ah]=header([s param.istr]); set(hh,'FontSize',10);
   % ss=[istr 10 ss];
  end

  [hh,ah]=header('fleft',ss);
  set(hh,'FontSize',8, 'VerticalAlign','bottom');

  if dy
    if ~isempty(fh) && isempty(get(hh,'UserData'))
       set(hh,'Pos', get(hh,'Pos')+[0 dy 0], 'UserData',1);
    end
  end

  if isfield(param,'sym') && ~isempty(param.sym)
     h=header('hleft');
     set(h,'String',[ get(h,'String') ' [' param.sym ']']);
  end

  if ~nargout, clear hh ah; end

end

