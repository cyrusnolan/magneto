classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model
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
        ellipse_a
        ellipse_b
        kp_M
        ki_M
        p0
    end
    
    methods
        function obj = parameters(model, gravitationalParameter, tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, safetySprintConst, safetyDampingRatio, trussLength, ...
                payloadRadius, massPayload, massTruss, initialAngularVelocity, tetherDelta, ...
                kp_M, zero_loc_M, p0)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here

            obj.model = model;
            obj.mu = gravitationalParameter;
            obj.k = tetherSpringConst;
            obj.b = tetherDampingRatio;
            obj.ktr = trussSpringConst;
            obj.btr = trussDampingRatio;
            obj.ks = safetySprintConst;
            obj.bs = safetyDampingRatio;
            obj.d0 = trussLength/2;
            obj.l0 = sqrt(payloadRadius^2 + obj.d0^2) + p0;
            obj.ls0 = sqrt(2)*payloadRadius;
            obj.mp = massPayload;
            obj.mtr = massTruss/6;
            obj.w = initialAngularVelocity;
            obj.ellipse_a = payloadRadius+tetherDelta;
            obj.ellipse_b = payloadRadius-tetherDelta;
            obj.kp_M = kp_M;
            obj.ki_M = kp_M * zero_loc_M;
            obj.p0 = p0;
        end
    end
end

