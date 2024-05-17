classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
        kp_amc
        ki_amc
        amc_on_time
        kp_tl
        kd_tl
        p0
    end
    
    methods
        function obj = parameters(gravitationalParameter, tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, safetySprintConst, safetyDampingRatio, trussLength, ...
                payloadRadius, massPayload, massTruss, initialAngularVelocity, tetherDelta, ...
                kp_amc, zero_loc_amc, p0, kd_tl, zero_loc_tl, amc_on_time)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here

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
            obj.mtr = massTruss/3;
            obj.w = initialAngularVelocity;
            obj.ellipse_a = payloadRadius+tetherDelta;
            obj.ellipse_b = payloadRadius-tetherDelta;
            obj.kp_amc = kp_amc;
            obj.ki_amc = kp_amc * zero_loc_amc;
            obj.amc_on_time = amc_on_time;
            obj.kd_tl = kd_tl;
            obj.kp_tl = kd_tl*zero_loc_tl;
            obj.p0 = p0;
        end
    end
end

