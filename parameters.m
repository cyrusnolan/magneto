classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k
        b
        mu
        mP
        mT
    end
    
    methods
        function obj = parameters(springConst, dampingRatio, ...
                gravitationalParameter, massPayload, massTruss)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.k = springConst;
            obj.b = dampingRatio;
            obj.mu = gravitationalParameter;
            obj.mP = massPayload;
            obj.mT = massTruss;
        end
    end
end

