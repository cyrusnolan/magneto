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
        mp
        mtr
        w
        le
        de
        re
    end
    
    methods
        function obj = parameters(gravitationalParameter, tetherSpringConst, tetherDampingRatio, trussSpringConst, ...
                trussDampingRatio, trussLength, payloadRadius, massPayload, numPayload, massTruss, ...
                initialAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.mu = gravitationalParameter;
            obj.k = tetherSpringConst;
            obj.b = tetherDampingRatio;
            obj.kc = trussSpringConst;
            obj.bc = trussDampingRatio;
            obj.d0 = trussLength/2;
            obj.l0 = sqrt(payloadRadius^2 + obj.d0^2);
            obj.mp = massPayload/numPayload;
            obj.mtr = massTruss/6;
            obj.w = initialAngularVelocity;
            sys_eqn = @(x) [2*obj.k*(x(1) - obj.l0)*sqrt(x(1)^2 - x(2)^2)/x(1) - obj.mp*obj.w^2*sqrt(x(1)^2 - x(2)^2);
                            2*obj.k*(x(1) - obj.l0)*x(2)/x(1) + obj.kc*(2*x(2) - 2*obj.d0) - obj.mtr*obj.w^2*x(2)];
            x0 = [obj.l0; obj.d0];
            solution = fsolve(sys_eqn, x0);
            obj.le = solution(1);
            obj.de = solution(2);
            obj.re = sqrt(obj.le^2 - obj.de^2);
        end
    end
end

