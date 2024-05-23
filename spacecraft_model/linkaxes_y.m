function linkaxes_y(ax)
% sets range(y) equal for all subplots
% chooses biggest range
    
    % gather data
    for i = length(ax):-1:1
        children = ax(i).Children;
        for j = 1:length(children)
            if isa(children(j), 'matlab.graphics.chart.primitive.Line')
                ydat = children(j).YData;
            end
        end
        ranges(i) = range(ydat);
        ymid(i) = min(ydat) + ranges(i)/2;
    end

    % set axes lims
    % axes limits are center of the range + maxrange/2
    maxrange = max(ranges);
    tick_range = maxrange*0.45;
    for i = 1:length(ax)
        low = ymid(i) - maxrange/2;
        high = ymid(i) + maxrange/2;
        ylim(ax(i), [low, high]);
        Y_ticks = [ymid(i) - tick_range, ymid(i), ymid(i) + tick_range];
        yticks(ax(i), Y_ticks);
        Y_tick_labels = {"- "+sprintf("%.2f", tick_range), sprintf("%0.3e", ymid(i)),"+ "+sprintf("%.2f", tick_range)};
        yticklabels(ax(i), Y_tick_labels);
    end
end