classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        l0
        BWc0
        BRac0
        BRbc0
        BRdc0
        BRec0
        BVac0
        BVbc0
        BVdc0
        BVec0
        simout
    end
    
    methods
        function obj = spacecraft(payloadRadius, ...
                trussLength, spacecraftAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.l0 = sqrt(payloadRadius^2+(trussLength/2)^2);
            b1 = [1;0;0]; b3 = [0;0;1]; % basis vectors
            obj.BWc0 = spacecraftAngularVelocity*b3;
            obj.BRac0 = payloadRadius*b1;
            obj.BRbc0 = -payloadRadius*b1;
            obj.BRdc0 = trussLength/2*b3;
            obj.BRec0 = -trussLength/2*b3;
            obj.BVac0 = cross(obj.BWc0,obj.BRac0);
            obj.BVbc0 = cross(obj.BWc0,obj.BRbc0);
            obj.BVdc0 = zeros(3,1);
            obj.BVec0 = zeros(3,1);
        end
        
        function obj = sim(obj, stopTime, gravity, params, orbit)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                stopTime {mustBePositive,mustBeReal,mustBeNumeric}
                gravity {mustBeA(gravity,"logical")}
                params {mustBeA(params,"parameters")}
                orbit {mustBeA(orbit,"orbit")}
            end
            NRco0 = orbit.NR;
            NVco0 = orbit.NV;
            Nh_orbit = orbit.Nh;
            Nb2 = -NRco0/norm(NRco0);
            Nb3 = Nh_orbit/norm(Nh_orbit);
            Nb1 = cross(Nb2,Nb3);
            NQB = [Nb1 Nb2 Nb3];
            NRac0 = NQB*obj.BRac0;
            NRbc0 = NQB*obj.BRbc0;
            NRdc0 = NQB*obj.BRdc0;
            NRec0 = NQB*obj.BRec0;
            NVac0 = NQB*obj.BVac0;
            NVbc0 = NQB*obj.BVbc0;
            NVdc0 = NQB*obj.BVdc0;
            NVec0 = NQB*obj.BVec0;
            NRao0 = NRac0+NRco0;
            NRbo0 = NRbc0+NRco0;
            NRdo0 = NRdc0+NRco0;
            NReo0 = NRec0+NRco0;
            NVao0 = NVac0+NVco0;
            NVbo0 = NVbc0+NVco0;
            NVdo0 = NVdc0+NVco0;
            NVeo0 = NVec0+NVco0;
            varInit = "StopTime = "+stopTime+"; " + ...
                "gravity = "+gravity+"; " + ...
                "k = "+params.k+"; " + ...
                "l0 = "+obj.l0+"; " + ...
                "b = "+params.b+"; " + ...
                "mu = "+params.mu+"; " + ...
                "ma = "+params.mP/2+"; " + ...
                "ma_inv = "+2/params.mP+"; " + ...
                "mb = "+params.mP/2+"; " + ...
                "mb_inv = "+2/params.mP+"; " + ...
                "md = "+params.mT/2+"; " + ...
                "md_inv = "+2/params.mT+"; " + ...
                "me = "+params.mT/2+"; " + ...
                "me_inv = "+2/params.mT+"; " + ...
                "NRao0 = "+mat2str(NRao0)+"; " + ...
                "NRbo0 = "+mat2str(NRbo0)+"; " + ...
                "NRdo0 = "+mat2str(NRdo0)+"; " + ...
                "NReo0 = "+mat2str(NReo0)+"; " + ...
                "NVao0 = "+mat2str(NVao0)+"; " + ...
                "NVbo0 = "+mat2str(NVbo0)+"; " + ...
                "NVdo0 = "+mat2str(NVdo0)+"; " + ...
                "NVeo0 = "+mat2str(NVeo0)+";";                     
            model = "dynamics";
            load_system(model)
            set_param(model, 'PreloadFcn', varInit)
            save_system(model)
            close_system(model)
            obj.simout = sim(model+".slx");
        end

        function plot(obj)
            figure(1)
            title("Distance between mass D and E (Truss Length)")
            t = obj.simout.tout;
            NRde = obj.simout.NRde;
            for i=length(NRde):-1:1
                rde(i) = norm(NRde(i,:));
            end
            plot(t,rde)
        end
    end
end
