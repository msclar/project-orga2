%% Load all files from folder and generate a plot of the surface for each one

function graficar()
    fnames = dir('output/*.out');
    numfids = length(fnames);
    for K = 1:numfids
        [X,delimiterOut]=importdata(sprintf('output/%s',fnames(K).name));
        
        figure('visible', 'off');
        surf(X);
        view(-74,40);
        zlim([305 400]);
        
        filename = str2num(fnames(K).name(2:end-4)) % esto usa que el archivo es T12434.out
        saveas(gca, sprintf('images/%d.png', filename));
    end
    
    % To plot all png saved as .gif, run:
    % convert -delay 10 -loop 0 $(ls -1v *.png) myimage.gif
