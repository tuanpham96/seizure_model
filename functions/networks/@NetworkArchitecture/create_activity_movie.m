function create_activity_movie(obj, vid_name, t_res, resp, t_in, inp, further_marked, color_opt)

if nargin < 8
    color_opt = 'EI';
end

% Set positions of axes and location of legend
switch obj.coord.type
    case {'RECT', 'LATTICE'}
        pos_ax1 = [0.05, 0.15, 0.5, 0.75];
        loc_ax1_lgnd = 'southoutside';
        
        pos_ax2 = [0.6, 0.65, 0.35, 0.25];
        pos_ax3 = [0.6, 0.35, 0.35, 0.25];
        pos_ax4 = [0.6, 0.05, 0.35, 0.25];
    case 'SPHERE'
        pos_ax1 = [0.1, 0.25, 0.8, 0.7];
        loc_ax1_lgnd = 'best';
        
        pos_ax2 = [0.05, 0.02, 0.25, 0.2];
        pos_ax3 = [0.35, 0.02, 0.25, 0.2];
        pos_ax4 = [0.65, 0.02, 0.25, 0.2];
end

% Create video writer
vid_write = VideoWriter(sprintf('%s.avi', vid_name));
open(vid_write);

% Title style
ttl_style = {'FontSize', 12, 'FontWeight', 'normal'};

% First frame
% Draw the population and the stimulated neurons
[x, y, mrkscl, gobj4lgdn] = obj.visualize_population('', further_marked);
gobj4lgdn(3).DisplayName = 'Received external input';

set(gcf, 'units', 'normalized', 'position', [0,0,1,1], 'color', 'w');
ax1 = gca; hold(ax1, 'on');
set(ax1, 'xcolor', 'none', 'ycolor', 'none', 'position', pos_ax1);
legend(ax1, gobj4lgdn, 'Location', loc_ax1_lgnd);

% Axes for plotting time series
general_style = {'unit', 'normalized', 'visible', 'off'};
ax2 = axes(general_style{:}, 'position', pos_ax2);
ax3 = axes(general_style{:}, 'position', pos_ax3);
ax4 = axes(general_style{:}, 'position', pos_ax4);
arrayfun(@(ax) hold(ax, 'on'), [ax1, ax2, ax3, ax4]);

% Plotting input time series
if size(inp, 1) > 1
    inp = normalize_minmax(inp(further_marked, :));
end
plot(ax2, t_in, inp, '-k', 'LineWidth', 1);

% Plotting EXC time series
resp_E = normalize_minmax(resp(obj.pop_ind.PE,:));
plot(ax3, t_res, resp_E, 'LineWidth', 1, 'Color', obj.colors.PE);

% Plotting INH time seriese
resp_I = normalize_minmax(resp(obj.pop_ind.PI,:));
plot(ax4, t_res, resp_I, 'LineWidth', 1, 'Color', obj.colors.PI);

% Overlay with no activity for 1st frame
sc_handle = scatter(ax1, x, y, 70*mrkscl, 'w', 'filled', 'o');
title(ax1, 't = 0 s', ttl_style{:});

% Vertical line representing time
time_lines = arrayfun(@(ax) plot(ax, [0,0], [0,1], ':k', ...
    'LineWidth', 0.7), [ax2, ax3, ax4]);

% Save to video
vid_frame = getframe(gcf);
writeVideo(vid_write, vid_frame);

% Discretize amplitude levels into discrete colors
num_lvl = 100;
sep_colors = false; 
if isa(color_opt, 'function_handle')
    cmap = color_opt(num_lvl);
elseif ischar(color_opt) 
    if strcmpi(color_opt, 'EI')
        cmap_PE = flipud(color2whitegradient(obj.colors.PE, num_lvl));
        cmap_PI = flipud(color2whitegradient(obj.colors.PI, num_lvl));
        PE = obj.pop_ind.PE;
        PI = obj.pop_ind.PI;
        sep_colors = true; 
    else
        try 
            color_opt = str2func(color_opt);
        catch
            error('The string "%s" is not an available colormap function', color_opt);
        end
        cmap = color_opt(num_lvl);
    end
end

if sep_colors
    res_amp = resp;
    res_rng = [min(resp(PE,:), [], 'all'), max(resp(PE,:), [], 'all')];
    res_amp(PE,:) = (resp(PE,:) - res_rng(1))/abs(diff(res_rng));
    res_rng = [min(resp(PI,:), [], 'all'), max(resp(PI,:), [], 'all')];
    res_amp(PI,:) = (resp(PI,:) - res_rng(1))/abs(diff(res_rng));
    res_lvl = round((num_lvl - 1) * res_amp) + 1;
else
    res_rng = [min(resp, [], 'all'), max(resp, [], 'all')];
    res_amp = (resp - res_rng(1))/abs(diff(res_rng));
    res_lvl = round((num_lvl - 1) * res_amp) + 1;
end

% Plot and write frame for each time point
for i = 1:length(t_res)
    ti = t_res(i);
    ri = res_lvl(:,i);
    
    % Overlay with the current amplitude
    delete(sc_handle);
    if i == 1 && sep_colors
        sc_handle = gobjects(2,1);
    end
    
    if sep_colors
        sc_handle(1) = scatter(ax1, x(PE), y(PE), 120*mrkscl, cmap_PE(ri(PE),:), 'o', 'filled');
        sc_handle(2) = scatter(ax1, x(PI), y(PI), 120*mrkscl, cmap_PI(ri(PI),:), 'o', 'filled');
    else 
        sc_handle = scatter(ax1, x, y, 70*mrkscl, cmap(ri,:), 'o', 'filled');
    end

    ti_sec = floor(ti); 
    ti_ms = 1000*(ti - ti_sec); 
    title(ax1, sprintf('t = %.0fs, %03.fms', ti_sec, ti_ms), ttl_style{:});
    legend(ax1, gobj4lgdn, 'Location', loc_ax1_lgnd);
    
    % Move vertical time line
    delete(time_lines);
    time_lines = arrayfun(@(ax) plot(ax, ti*[1,1], [0,1], ':k', ...
        'LineWidth', 0.7), [ax2, ax3, ax4]);
    
    % Save to video
    vid_frame = getframe(gcf);
    writeVideo(vid_write, vid_frame);
end

close(vid_write);
end