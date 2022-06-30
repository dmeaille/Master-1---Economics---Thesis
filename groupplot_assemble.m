function groupplot_assemble(leg,names,tex,var1,var2, var3, var4)
% This function aims at combining the irf graphs for different outputs 
% from dynare models., and display them in one figure. 
% One recquirement is that all the structure array given as input, from 
% oo_.irfs, have the same length, and the same order for variable names
% To use the latex option, include 
% names = M_.endo_names;
% tex = M_.endo_names_tex;
% at the end of your mod file, and use them while calling the function.
% To include a legende in the inputs, you can have it written in a cell
% array of string, like for instance :
% C = {'leg1', 'leg2', ...}


% % Setting the color 
color1 = [0, 0, 0];
color3 = [0 0.5 0.7410];
color2 = [0.2, 0.0, 0.8];
color4 = [0,0.4470, 0.7410];
col_line = [0.6, 0.6, 0.8];

%% rajouter la possibilité de dire : si tex, alors on va récupérer les tex par correspondance, from M:8;endo_names_tex

mat1 = cell2mat(struct2cell(var1));
titles = extractBefore(fieldnames(var1), '_');

ndsets = 3;
if nargin == 5
    mat3 = nan *  mat1;
    mat4 = nan * mat1;
    mat2 = cell2mat(struct2cell(var2));
    ndsets = 2;
elseif nargin == 6
    mat4 = nan * mat1;
    mat2 = cell2mat(struct2cell(var2));
    mat3 = cell2mat(struct2cell(var3));
    ndsets = 3;
elseif nargin == 7
    mat2 = cell2mat(struct2cell(var2));
    mat3 = cell2mat(struct2cell(var3));
    mat4 = cell2mat(struct2cell(var4));
    ndsets = 4;
elseif nargin >= 8
    error ('Combine takes maximum 7 arguments, with only 4 datasets maximum')
elseif nargin <= 4
    error('Combine takes minimum 5 arguments, with at least 2 datasets')
end 

% --------- To calibrate according to the size of the data ---------% 
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


figure('Renderer', 'painters', 'Position', [0 0 700 1000]);
set(gcf,'color','w');
for i = 1:nvars
    subplot(nrows,ncols,i);
    plot(mat1(i,:), '-', 'LineWidth', 1, 'Color', color1); hold on
    plot(mat2(i,:), '-', 'LineWidth', .75, 'Color', color2); hold on
    plot(mat3(i,:), '-.', 'LineWidth', 1, 'Color', color3); hold on
    plot(mat4(i,:), '--', 'LineWidth', 0.75, 'Color', color4); hold on
    line([0 columns(mat1)], [0 0], 'Color', col_line)
    xlim([1 size(mat1,2)]);
    title(strcat('$',tex(find(strcmp(titles(i), names))),'$'), 'interpreter', 'latex')
end
tightfig;
legend(leg, 'interpreter', 'latex')

end
