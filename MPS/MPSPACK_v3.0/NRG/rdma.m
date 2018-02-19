% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

% if isfield(param,'JJ')
%    banner('TEST rerun setup2CKondo'); % to save running cluster job ;)
%    setup2CKondo
% end

  if isset('Inrg') && isfield(Inrg,'Lambda') % isset('NRG')
     Lambda=Inrg.Lambda; NRG=FOUT{end};
  else
     if exist('user_fout','var') && ~isempty(user_fout)
     NRG=user_fout; else NRG='./NRG/NRG'; end
     
     fn=[NRG '_info.mat'];
     load2(fn,'Lambda','param','E0');
     if ~exist('FC','var') || ~exist('Z','var')
        q=load2(fn,'ops'); structexp(q.ops,'FC','Z');
     end
     if ~isempty(whos('-FILE',fn,'zflags')), load2(fn,'zflags'), end
     if ~isempty(whos('-FILE',fn,'cflags')), load2(fn,'cflags'), end
     if ~isempty(whos('-FILE',fn,'Gamma' )), load2(fn,'Gamma' ), end
  end

% if ~exist('Gamma','var'),  disp('Need Gamma for plot.'); return; end
  if ~exist('Lambda','var')
  disp('Need Lambda for rescaling.'); return; end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

  if exist('noDMA')~=1, noDMA=0; end
  if ~exist('plotflag','var'), plotflag=1; end

  if ~noDMA

  % options: 'calcOps', 'NRho', 50, 'hof', 'bul', 'fDM', 'symDMA'
  % odma={ 'partial' }; setopts(odma,'-PARTIAL?');
    odma=setopts('-','-partial?','-PARTIAL?');
    % Idma.A4 and Idma.A4f makes up ~97% of Idma structure!
    % Wb,Jul29,13

  % if exist('bul','var') && bul, odma{end+1}='bul'; end
  % if exist('hof','var') && hof, odma{end+1}='hof'; end
    if isset('dmaVersion')
         odma{end+1}=dmaVersion;
    else odma{end+1}='fDM'; end

    setopts(odma,'T?','NRho?','zflags?','cflags?','rhoNorm?','nlog?',...
      'mspec?','emin?','emax?',...
      '-calcOps?','-nostore?','-locRho?','-calcRho?','-keven?','-kodd?');

    if ~exist('Z0','var'),
    Z0=Z; end % H0 equals a Wilson site

    if 0
      if ~exist('op1','var'), op1=[]; end
      if ~exist('op2','var'), op2=FC; end
    else
      if ~exist('op1','var'), op1=[FC FN]; end
      if ~exist('op2','var'), op2=[FC FC]; end
    end

  % A0 is first MPS basis set!
  % ===================================================== %
    [om,a0,Idma] = fdmNRG_QS(NRG,op1,op2,Z0, odma{:});
  % ===================================================== %
  % dmNRG => fdmNRG_QS
  % NB! NRGWilsonQS/fdmNRG_QS also switches to col-major / LRs convention!
  % Wb,Aug11,10

    a=a0; a(find(isnan(a)))=0; a=sum(a);
    printf('\nSum raw spectral data: %s (%.3g)', vec2str(a), max(abs(a-1)))
    printf('See structure Idma for more info.\n');

  end

% -------------------------------------------------------------------- %
  if plotflag, dma_plot, end
% -------------------------------------------------------------------- %

