classdef FitGuiDataClass < handle

    properties
        x = [];
        dx = [];
        y = [];
        dy = [];
        
        Fit = struct('x_fit', [], 'y_fit', [], 'a', [], 'aerr', [], 'chisq', [], 'RChiSquare', [], 'x_fit_plot', [], 'y_fit_plot', []);
    end
    
    
    
end
