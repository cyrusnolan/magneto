classdef parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mu
        k
        b
        kc
        bc
        d0
        l0
        ma
        mb
        md
        me
        mf
        mg
        w
    end
    
    methods
        function obj = parameters(gravitationalParameter,tetherSpringConst, tetherDampingRatio, trussSpringConst, ...
                trussDampingRatio, trussLength, tetherLength, massPayload, massTruss, spacecraftAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.mu = gravitationalParameter;
            obj.k = tetherSpringConst;
            obj.b = tetherDampingRatio;
            obj.kc = trussSpringConst;
            obj.bc = trussDampingRatio;
            obj.d0 = trussLength/2;
            obj.l0 = tetherLength;
            obj.ma = massPayload/4;
            obj.mb = massPayload/4;
            obj.mf = massPayload/4;
            obj.mg = massPayload/4;
            obj.md = massTruss/6;
            obj.me = massTruss/6;
            obj.w = spacecraftAngularVelocity;
        end
    end
end

