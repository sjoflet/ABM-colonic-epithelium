function [choice] = myPlotColors(num)
colorMat = [ 0 0 0 ; ... % 1 black
    1 1 1 ; % 2 white 
    1 0 0 ; ... % 3 red
    0 1 0 ; ... % 4 bright green
    0 0 1 ; ... % 5 bright blue
    0 1 1 ; ... % 6 cyan
    1 0 1 ; ... % 7 magenta

    %% Pastels
    0.94 0.67 0.67 ; ... % 8 salmon
    0.76 0.86 0.61 ; ... % 9 light olive green
    0.47 0.67 0.19 ; ... % 10 dark olive green 
    0.31 0.81 0.81 ; ... % 10 pastel cyan
    0.85 0.67 0.92 ; ... % 11 pastel purple

    %% Off-standard
    0 0.4470 0.7410 ; ... % 12 cool blue
    0.6350 0.0780 0.1840 ;...  % 13 marroon
    0.8500 0.3250 0.0980 ; ... % 14 burnt orange
    0.9290 0.6940 0.1250 ; ... % 15 pollen yellow
    0.4940 0.1840 0.5560 ; ... % 16 purple
    0.4660 0.6740 0.1880 ; ... % 17 green
    0.3010 0.7450 0.9330 ; ... % 18 sky blue
    1 1 0 ] ; % 19 yellow

choice = [colorMat(num,:)] ;

end 

