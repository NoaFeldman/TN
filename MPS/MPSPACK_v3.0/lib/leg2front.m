function lh=leg2front()
% bring legend to the front of the axis set again
% Wb,Sep20,05

  ah=gca; lh=find_legend(ah);

  if ~isempty(lh), axes(lh); setax(ah); else
     lh=legend('toggle'); setax(lh); legend boxoff 
     if ~nargout, clear lh; end
  end

return

  m=0; ah=gca;
  for h=[ findall(gcf,'tag','legdetached'), findall(gcf,'tag','legend') ]'
     if isequal(ah,getuser(h,'parent')), axes(h); m=m+1; end
  end

end

% modified from matlab legend::find_legend() subroutine
% Wb,Jul29,11

function lh=find_legend(ah)

  lh=[]; if isempty(ah) || ~ishandle(ah), return; end

  lh=getappdata(ah,'LegendPeerHandle');

  if isempty(lh) || ~ishandle(lh) || ~isequal(get(lh,'Axes'),handle(ah))
    % reset the "LegendPeerHandle" appdata
    % rmappdata(ah,'LegendPeerHandle');
      lh=[];
  end

end

