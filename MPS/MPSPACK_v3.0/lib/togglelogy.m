function togglelogy(ah)
% Function - togglelogy
%
%    toggle scale of y-axis between log|y| and lin(y)
%    stores original lin data in UserData of axis set.
%
% Wb,Oct28,05

  if nargin==0, ah=gca; end
  makelog=get(ah,'YScale'); makelog=~isequal(makelog,'log');

  if makelog
      s=get(ah,'UserData'); if ~isstruct(s), clear s; end
      s.xlm = get(ah,'XLim');
      s.ylm = get(ah,'YLim');
      set(ah,'UserData',s);

      for h=findall(ah,'Type','Line')'
          l=get(h,'UserData'); if ~isstruct(l), clear l; end
          l.yd=get(h,'YData');
          set(h,'YData',abs(l.yd),'UserData',l);
      end

      axis tight % important to avoid warning on negative data (!???)
      set(ah,'YLimMode','auto','YScale','log'); axis tight
      xlim(s.xlm); yl=ylim; yl(2)=max(abs(s.ylm)); ylim(yl);
  else
      set(ah,'YScale','lin');
      s=get(ah,'UserData'); e=0;
      if all(isfield(s,{'ylm','ylm'}))
         set(ah,'XLim', s.xlm, 'YLim', s.ylm);

         for h=findall(ah,'Type','Line')'
             l=get(h,'UserData');
             if ~isfield(l,{'yd'}), e=e+1; continue; end
             set(h,'YData', l.yd);
         end
         if e
         wblog('ERR No or invalid UserData to return to lin-scale! (%d)',e); end
     end
  end

  if exist('xmark','file')
     xmark -fix
     ymark -fix
  end

end

