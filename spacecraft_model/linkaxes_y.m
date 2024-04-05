function linkaxes_y(ax)
% sets range(y) equal for all subplots
% chooses biggest range
% can handle case with multple y axes on a single subplot
    
    % gather ranges
    for i = length(ax):-1:1
        for j = length(ax(i).YAxis):-1:1
            % i = number of subplots
            % j = number of axes on subplot
            ranges(j, i) = range(ax(i).YAxis(j).Limits);
        end
    end

    % set axes lims
    for j = length(ax(i).YAxis):-1:1
        maxRange = max(ranges(j, :));
        for i = 1:length(ax)
            meanLimVal = mean(ax(i).YAxis(j).Limits);
            ax(i).YAxis(j).Limits = [meanLimVal-maxRange/2  meanLimVal+maxRange/2];
        end
    end
end

