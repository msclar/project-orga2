%% Load all files from folder and generate a plot of the surface for
%% each one

function graficar()
    fnames = dir('output/*.out');
    numfids = length(fnames);
    for K = 1:numfids
        fnames(K).name
        sprintf('output/%s',fnames(K).name)
        
        [X,delimiterOut]=importdata(sprintf('output/%s',fnames(K).name));
        
        f = figure('visible', 'on');
        surf(X);
        
        saveas(gca, sprintf('images/T%d.fig', K));
    end
