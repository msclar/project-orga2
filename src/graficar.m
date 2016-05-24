function graficar()


%% Load all files from folder
    fnames = dir('output/*.out');
    numfids = length(fnames);
    %vals = cell(1,numfids);
    for K = 1:numfids
        fnames(K).name
        sprintf('output/%s',fnames(K).name)
        
        [X,delimiterOut]=importdata(sprintf('output/%s',fnames(K).name));
        
        f = figure('visible', 'on');
        surf(X);
        %savefig(f, sprintf('images/T%d.fig', K));
        %close(f)

        %imwrite(surf(X), sprintf('images/T%d.png', K)); % 'images/T' num2str(K) '.png');
        %savefig(sprintf('images/T%d.fig', K)); % 'images/T' num2str(K) '.png');
        saveas(gca, sprintf('images/T%d.fig', K));
        %break;
    end

 
    %figure
    %surf(X);