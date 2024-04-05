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
        invIm
        rm
        w
    end
    
    methods
        function obj = parameters(model, gravitationalParameter, tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, safetySprintConst, safetyDampingRatio, trussLength, ...
                payloadRadius, massPayload, massTruss, Imotor, rmotor, initialAngularVelocity)
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
            obj.l0 = sqrt(payloadRadius^2 + obj.d0^2);
            obj.ls0 = sqrt(2)*payloadRadius;
            obj.mp = massPayload;
            obj.mtr = massTruss/6;
            obj.invIm = Imotor;
            obj.rm = rmotor;
            obj.w = initialAngularVelocity;
        end
    end
end

