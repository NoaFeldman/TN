function deltaVsW()
    % I have prepared data for longer systems, but they seem to divert form
    % CFT, so until I prepare mpre accurate simulations, I'll stick to
    % this length scale.
    Ls = 112:40:432;
    deltas = [-0.5 -0.3 -0.1 0 0.1 0.3 0.5 0.7 0.9 1];
    deltaStrings = {'0-5', '0-3', '0-1', '00', '01', '03', '05', '07', '09', '1'};
    for i = 1:length(deltas)
        hold off
        delta = deltas(i);
        spec = processRDMSpectrum(strcat('groundStateSpectrumDMRG_',...
            int2str(Ls(length(Ls))), '_', int2str(Ls(length(Ls))/2), '_', deltaStrings{i}));
        [as, aerrs, covs, chi2s, yfit] = fitnonlin(spec.alphas.', spec.alphas.', real(spec.salpha).', 0.01.*spec.alphas.', ...
            0.01.*real(spec.salpha).', 'sAlphaGSTheo', [3200 0.95 1], [1 K(delta)]);
        plot(spec.alphas.', real(spec.salpha).');
        hold on
        plot(spec.alphas.', yfit);
        title(strcat('$\Delta = ', num2str(delta), '$'), 'Interpreter', 'latex');
        X(i) = as(1);
        w(i) = as(2);
        pause(0.5);
    end
    hold off
    
    plot(deltas, X);
    
    plot(deltas, w);
    
    % X is not a constant and his is weird - I would expect X to equal L*c,
    % and I'm keeping L = 432 c onstant here.
    % I can try to fix X, and see where it gets me:
    fixedX = X(4); % X(delta = 0)
    for i = 1:length(deltas)
        hold off
        delta = deltas(i);
        spec = processRDMSpectrum(strcat('groundStateSpectrumDMRG_',...
            int2str(Ls(length(Ls))), '_', int2str(Ls(length(Ls))/2), '_', deltaStrings{i}));
        [as, aerrs, covs, chi2s, yfit] = fitnonlin(spec.alphas.', spec.alphas.', real(spec.salpha).', 0.01.*spec.alphas.', ...
            0.01.*real(spec.salpha).', 'sAlphaGSTheo_fixedX', [0.95 1], [1 K(delta) fixedX]);
        plot(spec.alphas.', real(spec.salpha).');
        hold on
        plot(spec.alphas.', yfit);
        title(strcat('$\Delta = ', num2str(delta), '$'), 'Interpreter', 'latex');
        w1(i) = as(1);
        pause(0.5);
    end
    plot(deltas, w1);
end

function K = K(jz)
    K = 1 ./ (2 * acos(-jz) / pi);
end