function [h, xscale, yvals] = skyline(varargin)
%SKYLINE Creates a skyline plot from single line plot data
%
% Syntax:
%  h = f.skyline(xdata, ydata);
%  [h, xscale, yvals] = f.skyline(ax, xdata, ydata, 'Name', value,...);
%
% Inputs:
%  ax - (def: gca) Can be given as first argument to specify target axes.
%  xdata - Line xdata vector
%  ydata - Line ydata vector (same number of points as xdata)
%  varargin - (Optional) 'Name',value pairs (same as built-in `plot`).
%
% Output:
%  h - Handle to skyline `line` graphics object
%  [xscale, yvals] - Values given to `plot` to produce `h`.
%
% See also: plot

if isa(varargin{1}, 'matlab.graphics.axis.Axes')
   ax = varargin{1};
   varargin(1) = [];
else
   ax = gca;
end
xdata = varargin{1};
ydata = varargin{2};
varargin(1:2) = [];
n = 2 * length(ydata);
yvals = zeros(1, n);
xscale = zeros(1, n);
yvals(1:2:n-1) = ydata;
yvals(2:2:n) = ydata;
xscale(1:2:n-1) = xdata;
xscale(2:2:n-2) = xdata(2:end);
xscale(end) = 2 * xdata(end) - xdata(end-1);
h = line(ax, xscale, yvals, ...
   'DisplayName', 'PSTH', ...
   'LineWidth', 2, ...
   'Tag', 'skyline', ...
   'Color', ax.ColorOrder(ax.ColorOrderIndex,:), ...
   varargin{:});
ax.ColorOrderIndex = mod(ax.ColorOrderIndex,7) + 1;
end