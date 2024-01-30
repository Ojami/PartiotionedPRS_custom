function ax = shade(ax, opts)
% shades input axis according to a desired shading length.
% Oveis Jamialahmadi, University of Gothenburg, March 2022.

arguments
    ax (1,1) {mustBeA(ax, 'matlab.graphics.axis.Axes')}
    opts.length (1,1) double {mustBeGreaterThan(opts.length, 1e-4)} = 1
    opts.step double = nan % to be used instead of spaced steps in 'length' (overrides 'length'). Must be a matrix of 2XN if dir is 'xy' and a vector otherwise.
    opts.FaceAlpha (1,1) double {mustBeInRange(opts.FaceAlpha, 1e-6, 1)} = 0.2
    opts.dir {mustBeTextScalar, mustBeMember(opts.dir, ["x", "y", "xy"])} = "y" % shades direction; x: vertical, y: horizontal, xy: both
end

xlim = ax.XLim; ylim = ax.YLim; % don't allow X/Y-lims to vary
if iscategorical(xlim)
    xlim = [0.5, numel(ax.XTick)+0.5];
end

if ~any(isnan(opts.step))
    if opts.dir == "xy"
        stpsx = opts.step(1, :);
        stpsy = opts.step(2, :);
    else
        stps = opts.step;
    end
else
    if opts.dir == "y"
        stps = ylim(1):opts.length:ylim(2);
    elseif opts.dir == "x"
        stps = xlim(1):opts.length:xlim(2);
    else
        stpsx = xlim(1):opts.length:xlim(2);
        stpsy = ylim(1):opts.length:ylim(2);
    end
end


if opts.dir ~= "xy" && numel(stps) < 3 % not enough points
    return
elseif opts.dir == "xy" && (numel(stpsx) < 3 || numel(stpsy) < 3)
    return
end

% keep only YData corresponding to shading areas
if opts.dir ~= "xy" 
    if rem(numel(stps), 2), stps(end) = []; end
    stps = reshape(stps, 2, [])';
elseif opts.dir == "xy"
    if rem(numel(stpsx), 2), stpsx(end) = []; end
    if rem(numel(stpsy), 2), stpsy(end) = []; end
    stpsx = reshape(stpsx, 2, [])';
    stpsy = reshape(stpsy, 2, [])';
end

if opts.dir == "x"
    for i = 1:size(stps, 1)
        patch(ax, 'YData', [ylim, flip(ylim)], ...
            'XData', repelem(stps(i, :), 2), ...
            'LineStyle', 'none', ...
            'FaceAlpha', opts.FaceAlpha, 'Tag', "shadefnc");
    end
elseif opts.dir == "y"
    for i = 1:size(stps, 1)
        patch(ax, 'XData', [xlim, flip(xlim)], ...
            'YData', repelem(stps(i, :), 2), ...
            'LineStyle', 'none', ...
            'FaceAlpha', opts.FaceAlpha, 'Tag', "shadefnc");
    end
else % xy
    for i = 1:size(stpsx, 1)
        patch(ax, 'YData', [ylim, flip(ylim)], ...
            'XData', repelem(stpsx(i, :), 2), ...
            'LineStyle', 'none', ...
            'FaceAlpha', opts.FaceAlpha, 'Tag', "shadefnc");
    end

    for i = 1:size(stpsy, 1)
        patch(ax, 'XData', [xlim, flip(xlim)], ...
            'YData', repelem(stpsy(i, :), 2), ...
            'LineStyle', 'none', ...
            'FaceAlpha', opts.FaceAlpha, 'Tag', "shadefnc");
    end
end

try ax.XLim = xlim; catch, end
try ax.YLim = ylim; catch, end % don't allow X/Y-lims to vary

pax = findobj(ax, 'Tag', "shadefnc");
pax_idx = ismember(ax.Children, pax);
ax.Children = [ax.Children(~pax_idx); ax.Children(pax_idx)]; % send patch to back so it doesn't cover other objects

end % END