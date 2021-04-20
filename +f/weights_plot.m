function [h, ax, fig] = weights_plot(pars, name)
%WEIGHTS_PLOT Creates a labeled weight matrix figure.
%
% Syntax:
%  h = f.weights_plot(pars);
%  [h, ax, fig] = f.weights_plot(pars, name);
%
% Inputs:

if ~isfield(pars, 'fig')
   if ~isfield(pars, 'ax')
      fig = figure(...
         'Name', 'Weights', ...
         'Position', [100 100 590 500], ...
         'Color', 'w',...
         'PaperOrientation', 'portrait',...
         'PaperType', 'usletter');
   else
      fig = pars.ax.Parent;
   end
else
   fig = pars.fig;
end

if ~isfield(pars, 'ax')
   ax = axes(fig,...
      'NextPlot', 'add', ...
      'XTick', 21:40:221, ...
      'YTick', 21:40:221,...
      'YTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'},...
      'XTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'});
else
   ax = pars.ax;
end

nColUnits = pars.n_excit + pars.n_inhib + pars.n_out;  % number of units in a column
nDispUnits = pars.n_excit + pars.n_inhib;  % number of units in col to display
nDisplay = 3 * (pars.n_excit + pars.n_inhib);  % Ingnore column output units.

a1Index = 1;
aOutIndex = a1Index + (pars.n_excit + pars.n_inhib);
b1Index = nColUnits + a1Index;
bOutIndex = b1Index + (pars.n_excit + pars.n_inhib);
c1Index = nColUnits + b1Index;
cOutIndex = c1Index + (pars.n_excit + pars.n_inhib);

wvals = zeros(nDisplay, nDisplay);
for iw = 1:pars.weights
   pre_unit = pars.weight_pre_unit(iw);
   post_unit = pars.weight_post_unit(iw);
   
   if pre_unit < b1Index
      ix = pre_unit;
      if pre_unit >= aOutIndex
         ix = 0;
      end
   elseif pre_unit < c1Index
      ix = pre_unit - pars.n_out;
      if pre_unit >= bOutIndex
         ix = 0;
      end
   else
      ix = pre_unit - 2 * pars.n_out;
      if pre_unit >= cOutIndex
         ix = 0;
      end
   end
   
   if post_unit < b1Index
      iy = post_unit;
      if post_unit >= aOutIndex
         iy = 0;
      end
   elseif post_unit < c1Index
      iy = post_unit - pars.n_out;
      if post_unit >= bOutIndex
         iy = 0;
      end
   else
      iy = post_unit - 2 * pars.n_out;
      if post_unit >= cOutIndex
         iy = 0;
      end
   end
   
   if (ix > 0) && (iy > 0)
      wvals(ix, iy) = pars.weight_strength(iw);
   end
end
wvals(wvals == 0) = NaN;
wvals(1, end) = pars.max_strength *  pars.PSP_factor;
wvals(end, 1) = -pars.max_strength *  pars.PSP_factor;
h = pcolor(ax, wvals/ pars.PSP_factor);

for itick = 41:40:241
   line(ax, [1 241], [itick itick], 'Color', 'k');
   line(ax, [itick itick], [1 241], 'Color', 'k');
end

ylabel(ax, 'Pre-synaptic Units');
xlabel(ax, 'Post-synaptic Units');
if nargin > 1
   title(ax, name, 'Color', 'k', 'FontSize', 16, ...
      'FontWeight', 'bold');
else
   title(ax, 'Weight Matrix', 'Color', 'k', 'FontSize', 16, ...
      'FontWeight', 'bold');
end
colormap jet;
shading flat;
caxis(ax, [-500 500]);
colorbar(ax);
drawnow
end