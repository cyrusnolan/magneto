classdef orbit
    %ORBIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NR
        NV
        Nh
        n
    end
    
    methods
        function obj = orbit(gravitationalParameter, semiMajorAxis, ...
                eccentricity, inclination, argumentPeriapsis, ...
                lonAscendingNode, trueAnomaly)
            %ORBIT Construct an instance of this class
            %   Detailed explanation goes here
            p = [1;0;0]; q = [0;1;0]; % perifocal basis vectors
            mu = gravitationalParameter;
            a = semiMajorAxis;
            e = eccentricity;
            i = inclination;
            w = argumentPeriapsis;
            W = lonAscendingNode;
            v = trueAnomaly;
            slr = a*(1-e^2);
            h = sqrt(mu*slr);
            PR = slr/(1+e*cos(v))*(cos(v)*p+sin(v)*q);
            PV = mu/h*(-sin(v)*p+(e+cos(v))*q);
            Ph = cross(PR,PV);
            NQF = [cos(W) -sin(W) 0; sin(W) cos(W) 0; 0 0 1];
            FQG = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
            GQP = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];
            NQP = NQF*FQG*GQP;
            obj.NR = NQP*PR;
            obj.NV = NQP*PV;
            obj.Nh = NQP*Ph;
            obj.n = sqrt(mu/a^3);
        end
    end
end

