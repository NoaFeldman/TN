function sFullVST(dirname) 
    specResults = arrangeDMRGResults(dirName, 0.5);
    fig = figure('visible','off');
    scatter(specResults.t, specResults.sFull);
    xticks(specResults.t(0:10:end))
    grid on
    saveas(fig,[dirname '/sFullVST'],'png');
end