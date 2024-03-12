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
        function obj = parameters(config, gravitationalParameter, tetherSpringConst, tetherDampingRatio, ...
                trussSpringConst, trussDampingRatio, trussLength, payloadRadius, massPayload, massTruss, ...
                initialAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                config {mustBeMember(config,['v1','octahedron'])}
                gravitationalParameter
                tetherSpringConst
                tetherDampingRatio
                trussSpringConst
                trussDampingRatio
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
            obj.d0 = trussLength/2;
            obj.l0 = sqrt(payloadRadius^2 + obj.d0^2);
            obj.mp = massPayload;
            obj.mtr = massTruss/6;
            obj.w = initialAngularVelocity;
            obj = setEquilibrium(obj);
        end

        function obj = setEquilibrium(obj)
            switch obj.config
                case "v1"
                    numPayload = 2;
                case "octahedron"
                    numPayload = 4;
            end
            sys_eqn = @(x) [2*obj.k*(x(1) - obj.l0)*sqrt(x(1)^2 - x(2)^2)/x(1) - obj.mp*obj.w^2*sqrt(x(1)^2 - x(2)^2);
                            numPayload*obj.k*(x(1) - obj.l0)*x(2)/x(1) + obj.ktr*(2*x(2) - 2*obj.d0)];
            x0 = [obj.l0; obj.d0];
            solution = fsolve(sys_eqn, x0);
            obj.le = solution(1);
            obj.de = solution(2);
            obj.re = sqrt(obj.le^2 - obj.de^2);
        end
    end
end

