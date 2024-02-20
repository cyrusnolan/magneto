classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k
        b
        kc
        bc
        mu
        mP
        mT
    end
    
    methods
        function obj = parameters(tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, ...
                gravitationalParameter, massPayload, massTruss)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.k = tetherSpringConst;
            obj.b = tetherDampingRatio;
            obj.kc = trussSpringConst;
            obj.bc = trussDampingRatio;
            obj.mu = gravitationalParameter;
            obj.mP = massPayload;
            obj.mT = massTruss;
        end
    end
end

