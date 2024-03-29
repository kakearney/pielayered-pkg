function h = pielayered(x, lbl, total)
%PIELAYERED Plot a subdivided pie chart
%
% h = pielayered(x, lbl, total)
%
% Input variables:
%
%   x:      vector, length nx, data for pie chart.  Always treated the same
%           even if sum(x) <= 1 (unlike pie) 
%
%   lbl:    vector cell array, each cell holds a cell array of strings
%           associated with each data point.  Each string indicates a
%           category associated with the data point, from top level to most
%           local level    
%
%   total:  number indicating a full pie.  If sum(x) < total, then the pie
%           will be less than a full circle.  Default to sum(x) if not
%           included.



% Note: labels w/ slashes won't work right now

if ~isvector(x)
    error('x should be vector');
end
x = x(:);
nx = length(x);

if ~iscell(lbl) || ~isvector(lbl) || ~isequal(nx, length(lbl))
    error('lbl should be vector cell array same size as x');
end
lbl = lbl(:);

if ~all(cellfun(@ischar, cat(2, lbl{:})))
    error('All labels must be strings');
end
if nargin < 3
    total = sum(x);
end

if sum(x) > total
    error('Total must be <= sum(x)');
end

for ii = 1:nx
    if ~iscell(lbl{ii})
        lbl{ii} = {lbl{ii}};
    end
end

maxlev = max(cellfun(@length, lbl));

tmp = cellfun(@(x) fullfile(x{:}), lbl, 'uni', 0);
[srt, isrt] = sort(tmp);
lbl = lbl(isrt);
x = x(isrt);

xtot = zeros(nx,maxlev);
for ilev = 1:maxlev
    
    lbltot(:,ilev) = cellfun(@(a) fullfile(a{1:min(length(a),ilev)}), lbl,'uni', 0);
   
    [unq, ui, uj] = unique(lbltot(:,ilev));
    xtot(ui,ilev) = splitapply(@sum, x, uj);
    % [blah, xtot(ui,ilev)] = consolidator(uj, x, @sum);
  
end

xtot = xtot./total;

% Calculate coordinates

maxpt = 100;
r = 0:maxlev;

plotdata = cell(0,7);
for ilev = 1:maxlev
    
    theta = [0; cumsum(xtot(:,ilev))*2*pi];
    for is = 1:nx
        if ilev == 1 || ~strcmp(lbltot{is,ilev}, lbltot{is,ilev-1})
            
            dth = theta(is+1) - theta(is);
            
            if ~(dth == 0)
                npt = max(2, ceil(dth/(2*pi/maxpt)));
                th = linspace(theta(is), theta(is+1), npt);

                xouter = r(end) .* cos(th);
                youter = r(end) .* sin(th);
                xinner = r(ilev) .* cos(th);
                yinner = r(ilev) .* sin(th);
                x = [xouter fliplr(xinner) xouter(1)];
                y = [youter fliplr(yinner) youter(1)];

                xc = (r(ilev)+r(ilev+1))./2 .* cos((theta(is+1) + theta(is))/2);
                yc = (r(ilev)+r(ilev+1))./2 .* sin((theta(is+1) + theta(is))/2);

                lblparts = regexp(lbltot{is,ilev}, filesep, 'split');
                plotdata = [plotdata; {is, ilev, x, y, xc, yc, lblparts{end}}];
            end
            
        end
    end
    
    
end
    
% Plot

for ip = 1:size(plotdata)
    h.p(ip) = patch(plotdata{ip,3}, plotdata{ip,4}, 'w');
    h.t(ip) = text(plotdata{ip,5}, plotdata{ip,6}, plotdata{ip,7}, 'horiz', 'center');
end    
uistack(h.t, 'top');
h.lev = plotdata(:,2);
h.slice = plotdata(:,1);
set(gca, 'visible', 'off', 'dataaspectratio', [1 1 1]);



