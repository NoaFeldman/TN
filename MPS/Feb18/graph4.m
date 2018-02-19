function data = graph4()
    filenames = ["spectrumN2000D0_1024.mat", "spectrumN2000D0.5_1024.mat", "spectrumN2000D0.9_1024.mat", "spectrumN1500D1_1024.mat"];
    Delta = [0, -0.5, -0.9, 1];
    for f = 1 : length(filenames)
        s = load(filenames(f));
        marr = s.spectrum.keys();
        sums = zeros(length(marr));
        nums = zeros(length(marr));
        totalSum = 0;
        totalNum = 0;
        for m = 1:length(marr)
            vals = s.spectrum(marr{m});
            % keep only eigenvalues > 1e-9
            for i = 1:length(vals)
                if vals(i) < 1e-9
                    vals = vals(1 : i-1);
                    break;
                end
                sums(m) = sums(m) + vals(i);
            end
            nums(m) = length(vals);
            if (str2num(marr{m}) == 0)
                totalSum = totalSum + sums(m);
                totalNum = totalNum + nums(m);
            else 
                totalSum = totalSum + 2 * sums(m);
                totalNum = totalNum + 2 * nums(m);               
            end
        end
        weights = containers.Map();
        statesNums = containers.Map();
        for m = 1:length(marr)
            sz = str2num(marr{m}) / 2;
            if (nums(m) ~= 0)
                weights(num2str(sz)) = sums(m) / totalSum;
                statesNums(num2str(sz)) = nums(m);
            end
        end
        data.(strcat('delta', num2str(f), 'Weights')) = weights;
        data.(strcat('delta', num2str(f), 'StatesNums')) = statesNums;        
        data.(strcat('delta', num2str(f), 'N')) = totalNum;
        data.(strcat('delta', num2str(f))) = struct;
        
        disp(['Delta = ' num2str(Delta(f)) ', N_lambda = ' int2str(totalNum)]);
        disp(['m' char(9) 'N' char(9) 'weight']);
        for k = keys(statesNums)
            m = k{1};
            sz = str2num(m);
            disp([num2str(sz) char(9) num2str(statesNums(m)) ...
                char(9) num2str(weights(m))]);
        end
    end
end