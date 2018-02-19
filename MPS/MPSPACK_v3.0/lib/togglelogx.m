function togglelogx(ah)
% Function - togglelogx
%
%    toggle scale of x-axis between log|x| and lin(x)
%    stores original lin data in UserData of axis set.
%
% Wb,Oct28,05

  if nargin==0, ah=gca; end
  makelog=get(ah,'XScale'); makelog=~isequal(makelog,'log');

  if makelog
      s=get(ah,'UserData'); if ~isstruct(s), clear s; end
      s.xlm = get(ah,'XLim');
      s.ylm = get(ah,'YLim'); set(ah,'UserData',s);

      for h=findall(ah,'Type','Line')'
          l=get(h,'UserData'); if ~isstruct(l), clear l; end
          l.xd=get(h,'XData');
          set(h,'XData', abs(l.xd), 'UserData', l);
      end

    % axis tight % important to avoid warning on negative data (!???)
      set(ah,'XLimMode','auto','XScale','log'); axis tight
      xl=xlim; xl(2)=max(abs(s.xlm)); xlim(xl); ylim(s.ylm);
  else
      set(ah,'XScale','lin');
      s=get(ah,'UserData'); e=0;
      if all(isfield(s,{'xlm','ylm'}))
         set(ah,'XLim', s.xlm, 'YLim', s.ylm);

         for h=findall(ah,'Type','Line')'
             l=get(h,'UserData');
             if ~isfield(l,{'xd'}), e=e+1; continue; end
             set(h,'XData', l.xd);
         end
         if e
         wblog('ERR No or invalid UserData to return to lin-scale!'); end
     end
  end

  if exist('xmark','file')
     xmark -fix
     ymark -fix
  end

end

