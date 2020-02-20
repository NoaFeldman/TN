colors = zeros(3, 4);
colors(:, 1) = [0 0 0.5];
colors(:, 2) = [0 0.7 0.5];
colors(:, 3) = [0 0.7 0];
colors(:, 4) = [0.7 0.7 0];
scatter(t, rQ0, '.', 'markeredgecolor', colors(:, 1));
hold on
scatter(t, rQ1, '.', 'markeredgecolor', colors(:, 2));
scatter(t, rQ2, '.', 'markeredgecolor', colors(:, 3));
scatter(t, rQ3, '.', 'markeredgecolor', colors(:, 4));
xlabel('$t[\hbar/J]$', 'Interpreter', 'latex')
ylabel(strcat('$R_', int2str(n), '(Q)/R_', int2str(n), '$'), 'Interpreter', 'latex');
title(strcat('tDMRG results, $\Delta = ', num2str(Delta), '$'), 'Interpreter', 'latex');
legend({'$\Delta Q = 0$', '$\Delta Q = 1$', '$\Delta Q = 2$', '$\Delta Q = 3$'}, 'Interpreter', 'latex', 'FontSize', 14);
box on
savefig(strcat('~/Pictures/neg3_num_delta', num2str(10*Delta), 'l', int2str(l)))
saveas(gcf, strcat('~/Pictures/neg3_num_delta', num2str(10*Delta), 'l', int2str(l), '.png'))
hold off
plot(t(1:10:length(t)), c0, 'color', colors(:, 1))
hold on
plot(t(1:10:length(t)), c1, 'color', colors(:, 2))
plot(t(1:10:length(t)), c2, 'color', colors(:, 3))
plot(t(1:10:length(t)), c3, 'color', colors(:, 4))
xlabel('$t[\hbar/J]$', 'Interpreter', 'latex')
ylabel(strcat('$R_', int2str(n), '(Q)/R_', int2str(n), '$'), 'Interpreter', 'latex');
title(strcat('CFT Prediction'), 'Interpreter', 'latex');
legend({'$\Delta Q = 0$', '$\Delta Q = 1$', '$\Delta Q = 2$', '$\Delta Q = 3$'}, 'Interpreter', 'latex', 'FontSize', 14);
ylim([0 0.8])
box on
savefig(strcat('~/Pictures/neg3_theo_delta', num2str(10*Delta), 'l', int2str(l)))
saveas(gcf, strcat('~/Pictures/neg3_theo_delta', num2str(10*Delta), 'l', int2str(l), '.png'))
