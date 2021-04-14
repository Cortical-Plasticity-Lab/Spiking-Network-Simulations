function spikenet_wfig(usetitle, filename)
% Creates a weight matrix figure for a spikenet.mat parameter file that
% the user selects in the open file dialog box.

datapath = 'c:\data\spikenet\spikenet50';  % path to store networks
old_dir = cd(datapath);

if nargin < 2
    fname = '';
else
    fname = filename;
end

% Check if file is a .mat file

if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.mat') == 0)
    % Request file from user
    [uifname, uipath] = uigetfile('*.mat', 'Select Spikenet file');
    if (uifname(1) == 0) || (uipath(1) == 0)
        return
    end
    fname = fullfile(uipath, uifname);

    % Check for obvious errors.
    if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.mat') == 0)
        disp(['Error: ' fname ' is not a Spikenet file']);
        return;
    end
end

% Load spikenet parameters.

try
    index = strfind(fname, '\');
    if isempty(index)
        fpath = '';
        fprefix = fname(1:end-4);
    else
        fpath = fname(1:index(end));
        fprefix = fname(index(end)+1:end-4);
    end
    loadstruct = load(fname);
    p = loadstruct.p;
catch
    disp(['Error: ' fname ' is not a Spikenet  file']);
    return;
end

% Display Weight Matrix
wdfig = figure;
set(wdfig, 'Position', [100 100 420 400]);

p.n_excit = 40;  % Number of units in each excitatory group
p.n_inhib = 40;  % Number of units in each inhibitory group.
p.n_out = 40;    % Number of final output units for each Column
nColUnits = p.n_excit + p.n_inhib + p.n_out;  % number of units in a column
nDispUnits = p.n_excit + p.n_inhib;  % number of units in col to display
nDisplay = 3 * (p.n_excit + p.n_inhib);  % Ingnore column output units.

a1Index = 1;
aOutIndex = a1Index + (p.n_excit + p.n_inhib);
b1Index = nColUnits + a1Index;
bOutIndex = b1Index + (p.n_excit + p.n_inhib);
c1Index = nColUnits + b1Index;
cOutIndex = c1Index + (p.n_excit + p.n_inhib);

wvals = zeros(nDisplay, nDisplay);
for iw = 1:p.weights
    pre_unit = p.weight_pre_unit(iw);
    post_unit = p.weight_post_unit(iw);
    
    if pre_unit < b1Index
        ix = pre_unit;
        if pre_unit >= aOutIndex
            ix = 0;
        end
    elseif pre_unit < c1Index
        ix = pre_unit - p.n_out;
        if pre_unit >= bOutIndex
            ix = 0;
        end       
    else
        ix = pre_unit - 2 * p.n_out;
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
        iy = post_unit - p.n_out;
        if post_unit >= bOutIndex
            iy = 0;
        end       
    else
        iy = post_unit - 2 * p.n_out;
        if post_unit >= cOutIndex
            iy = 0;
        end               
    end

    if (ix > 0) && (iy > 0)
        wvals(ix, iy) = p.weight_strength(iw);
    end
end
wvals(wvals == 0) = NaN;
%wvals(1, end) = p.max_strength *  p.psp_factor;
%wvals(end, 1) = -p.max_strength *  p.psp_factor;
h = pcolor(wvals/ p.psp_factor);
caxis([-500 500]);

for itick = 41:40:241
    line([1 241], [itick itick], 'Color', 'k');
    line([itick itick], [1 241], 'Color', 'k');
end

set(gca, 'fontsize', 16);
set(gca, 'fontweight', 'bold');
set(gca, 'XTick', 21:40:221);
set(gca, 'YTick', 21:40:221);
set(gca, 'YTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'});
set(gca, 'XTickLabel', {'Ae', 'Ai', 'Be', 'Bi', 'Ce', 'Ci'});

ylabel('Pre-synaptic Units');
xlabel('Post-synaptic Units');
if nargin > 0
    title(usetitle)
else
    title('Weight Matrix');
end
colormap jet;
shading flat;
drawnow

cd(old_dir);