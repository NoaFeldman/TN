h = findobj(gca);
legend([h(6) h(2)], {'tDMRG', 'CFT'}, 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$t[\hbar/J]$', 'Interpreter', 'latex')
ylabel(strcat('$R_', int2str(n), '(Q)/R_', int2str(n), '$'), 'Interpreter', 'latex');
title(strcat('Charge resolved R\''enyi negativity, $\Delta = ', num2str(Delta), '$'), 'Interpreter', 'latex');
annotation('textbox',[.2 .2 .2 .6],'String', '$\Delta Q = 0$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
annotation('textbox',[.2 .2 .2 .4],'String', '$\Delta Q = 1$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
annotation('textbox',[.2 .2 .2 .2],'String', '$\Delta Q = 2$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
annotation('textbox',[.2 .2 .2 .0],'String', '$\Delta Q = 3$','FitBoxToText','on', 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 14);
box on
savefig(strcat('~/Pictures/neg_res_delta', num2str(10*Delta), 'l', int2str(l)));
saveas(gcf, strcat('~/Pictures/neg_res_delta', num2str(10*Delta), 'l', int2str(l), '.png'));