classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        config
        mu
        k
        b
        ktr
        btr
        ks
        bs
        d0
        l0
        ls0
        mp
        mtr
        w
    end
    
    methods
        function obj = parameters(config, gravitationalParameter, tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, safetySprintConst, safetyDampingRatio, trussLength, ...
                payloadRadius, massPayload, massTruss, initialAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                config
                gravitationalParameter
                tetherSpringConst
                tetherDampingRatio
                trussSpringConst
                trussDampingRatio
                safetySprintConst
                safetyDampingRatio
                trussLength
                payloadRadius
                massPayload
                massTruss
                initialAngularVelocity
            end
            obj.config = config;
            obj.mu = gravitationalParameter;
            obj.k = tetherSpringConst;
            obj.b = tetherDampingRatio;
            obj.ktr = trussSpringConst;
            obj.btr = trussDampingRatio;
            obj.ks = safetySprintConst;
            obj.bs = safetyDampingRatio;
            obj.d0 = trussLength/2;
            obj.l0 = sqrt(payloadRadius^2 + obj.d0^2);
            obj.ls0 = sqrt(2)*payloadRadius;
            obj.mp = massPayload;
            obj.mtr = massTruss/6;
            obj.w = initialAngularVelocity;
        end
    end
end

