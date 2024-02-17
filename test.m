

xorange = max(xco) - min(xco);
yorange = max(yco) - min(yco);
zorange = max(zco) - min(zco);
xomin = min(xco) - 0.1*xorange;
xomax = max(xco) + 0.1*xorange;
yomin = min(yco) - 0.1*yorange;
yomax = max(yco) + 0.1*yorange;
zomin = min(zco) - 0.1*zorange;
zomax = max(zco) + 0.1*zorange;

[oranges,idx] = sort([xorange yorange zorange]);
if idx(1) == 1
    xlim([min(xco)-0.1*oranges(2), max(xco)+0.1*oranges(2)])
    ylim([yomin yomax])
    zlim([zomin zomax])
elseif idx(1) == 2
    xlim([xomin xomax])
    ylim([min(yco)-0.1*oranges(2), max(yco)+0.1*oranges(2)])
    zlim([zomin zomax])
elseif idx(1) == 3
    xlim([xomin xomax])
    ylim([yomin yomax])
    zlim([min(zco)-0.1*oranges(2), max(zco)+0.1*oranges(2)])
end