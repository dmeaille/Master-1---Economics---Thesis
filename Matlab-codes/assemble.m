function combine_irfs_dynare(var, names, tex)
% This function aims at combining the irf graphs for different outputs 
% from dynare models. 
% One recquirement is that all the structure array given as input, from 
% oo_.irfs, have the same length, and the same order for variable names
% To use the latex option, include 
% names = M_.endo_names;
% tex = M_.endo_names_tex;
% at the end of your mod file, and use them while calling the function.
% To include a legende in the inputs, you can have it written in a cell
% array of string, that is using :
% C = {}

% % Setting the color 

color1 = [0.2, 0.1, 0.6];
color2 = [0.6, 0.6, 0.8];

if nargin == 1
    assignin('base','names',M_.endo_names);
    assignin('base', 'tex', M_.endo_names_tex)
end

%% rajouter la possibilité de dire : si tex, alors on va récupérer les tex par correspondance, from M:8;endo_names_tex

mat1 = cell2mat(struct2cell(var));
titles = extractBefore(fieldnames(var), '_');

% --------- Pour calibrer selon la taille des données ---------% 
nvars = size(titles,1);
if nvars <= 1 
    nrows= 3;
    ncols = 3;
elseif nvars <= 12
    nrows = 4;
    ncols = 3;
elseif nvars <= 15
    nrows = 5;
    ncols = 3;
elseif nvars <= 18
    nrows = 6;
    ncols = 3;
elseif nvars <= 21
    nrows = 7;
    ncols = 3;
elseif nvars <= 24
    nrows = 6;
    ncols = 4;
elseif nvars <= 28
    nrows = 7;
    ncols = 4;
else
    error('Too many variables')
end

% figure('NumberTitle', 'off', 'Name', 'This is the figure title');


figure('Renderer', 'painters', 'Position', [0 0 700 1000]);
set(gcf,'color','w');
for i = 1:nvars
    subplot(nrows,ncols,i);
    plot(mat1(i,:), '-', 'LineWidth', 1, 'Color', color1); hold on
    line([0 columns(mat1)], [0 0], 'Color', color2, "linewidth", 1.25);
    xlim([1 size(mat1,2)]);
    title(strcat('$',tex(find(strcmp(titles(i), names))),'$'), 'interpreter', 'latex')
end
tightfig;
end
