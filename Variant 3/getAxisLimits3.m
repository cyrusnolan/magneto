function limits = getAxisLimits3(data)
%GETAXISLIMITS Summary of this function goes here
%   Detailed explanation goes here
    arguments
        data (:,3) {mustBeNumeric, mustBeReal}
    end
    
    xmin = min(data(:,1));
    xmax = max(data(:,1));
    ymin = min(data(:,2));
    ymax = max(data(:,2));
    zmin = min(data(:,3));
    zmax = max(data(:,3));
    limits = [xmin xmax ymin ymax zmin zmax];

end

