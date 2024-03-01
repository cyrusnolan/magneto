classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        d0
        l0
        de
        re
        ma
        mb
        md
        me
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
        function obj = spacecraft(payloadRadius, massPayload, ...
                trussLength, massTruss, spacecraftAngularVelocity)
            %PARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            obj.d0 = trussLength/2;
            obj.l0 = sqrt(payloadRadius^2+obj.d0^2);
            obj.ma = massPayload/2;
            obj.mb = massPayload/2;
            obj.md = massTruss/6;
            obj.me = massTruss/6;
            b1 = [1;0;0]; b3 = [0;0;1]; % basis vectors
            obj.BWc0 = spacecraftAngularVelocity*b3;
            obj.BRac0 = payloadRadius*b1;
            obj.BRbc0 = -payloadRadius*b1;
            obj.BRdc0 = obj.d0*b3;
            obj.BRec0 = -obj.d0*b3;
            obj.BVac0 = cross(obj.BWc0,obj.BRac0);
            obj.BVbc0 = cross(obj.BWc0,obj.BRbc0);
            obj.BVdc0 = zeros(3,1);
            obj.BVec0 = zeros(3,1);
        end

        % function obj = getEquilibrium(obj, params)
        %     syms deq req
        %     we = norm(obj.BWc0);
        %     k = params.k;
        %     b = params.b;
        %     kc = params.kc;
        %     bc = params.bc;
        %     ma = params.ma;
        %     md = params.md;
        %     eqn1 = 2*k*(sqrt(req^2+deq^2)-obj.l0)*req/obj.l0 == ma*we^2*req;
        %     eqn2 = 2*k*(sqrt(req^2+deq^2)-obj.l0)*deq/obj.l0 == md*we^2*deq+2*kc*(deq-obj.d0);
        % 
        % end

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
                "b = "+params.b+"; " + ...
                "kc = "+params.kc+"; " + ...
                "bc = "+params.bc+"; " + ...
                "d0 = "+obj.d0+"; " + ...
                "l0 = "+obj.l0+"; " + ...
                "mu = "+params.mu+"; " + ...
                "ma = "+obj.ma+"; " + ...
                "ma_inv = "+1/obj.ma+"; " + ...
                "mb = "+obj.mb+"; " + ...
                "mb_inv = "+1/obj.mb+"; " + ...
                "md = "+obj.md+"; " + ...
                "md_inv = "+1/obj.md+"; " + ...
                "me = "+obj.me+"; " + ...
                "me_inv = "+1/obj.me+"; " + ...
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
            % find extrema
            NRao = obj.simout.NRao;
            NRbo = obj.simout.NRbo;
            for i=length(NRao):-1:1
                rao = norm(NRao(i,:));
                rbo = norm(NRbo(i,:));
                diff(i) = rao-rbo;
            end
            diffMinCond = [false, (diff(2:end-1) > diff(1:end-2)) & (diff(2:end-1) > diff(3:end)), false];
            diffMaxCond = [false, (diff(2:end-1) < diff(1:end-2)) & (diff(2:end-1) < diff(3:end)), false];
            diff0Cond = [false, (abs(diff(2:end-1)) < abs(diff(1:end-2))) & (abs(diff(2:end-1)) < abs(diff(3:end))), false];
            idxMinDiff = find(diffMinCond);
            idxMaxDiff = find(diffMaxCond);
            idx0Diff = find(diff0Cond);

            figure(1) % truss length (axes=axL) and tether length (axes=axR) vs time
            axL = subplot(1,2,1);
            hold on
            t = obj.simout.tout;
            rde = obj.simout.rde;
            plot(axL,t,rde)
            plot(axL,t(idxMinDiff),rde(idxMinDiff),'ro',LineWidth=2)
            plot(axL,t(idxMaxDiff),rde(idxMaxDiff),'bo',LineWidth=2)
            plot(axL,t(idx0Diff),rde(idx0Diff),'ko',LineWidth=2)
            title("Distance between mass D and E (Truss Length) vs time")
            xlabel("t (s)")
            ylabel("truss length (m)")
            legend('','local vertical alignment, payload a low','local vertical alignment, payload a high','local horizonal alignment')
            axR = subplot(1,2,2);
            hold on
            plot(axR,t,obj.simout.rad,'b')
            plot(axR,t,obj.simout.rae,'m-.')
            plot(axR,t,obj.simout.rbd,'k')
            plot(axR,t,obj.simout.rbe,'r-.')
            plot(axR,t(idxMinDiff),obj.simout.rad(idxMinDiff),'ro',LineWidth=2)
            plot(axR,t(idxMaxDiff),obj.simout.rad(idxMaxDiff),'bo',LineWidth=2)
            plot(axR,t(idx0Diff),obj.simout.rad(idx0Diff),'ko',LineWidth=2)
            plot(axR,t(idxMinDiff),obj.simout.rbd(idxMinDiff),'ro',LineWidth=2)
            plot(axR,t(idxMaxDiff),obj.simout.rbd(idxMaxDiff),'bo',LineWidth=2)
            plot(axR,t(idx0Diff),obj.simout.rbd(idx0Diff),'ko',LineWidth=2)
            title("Tether length vs time")
            xlabel("t (s)")
            ylabel("tether length (m)")
            legend('rad','rae','rbd','rbe','local vertical alignment, payload a low','local vertical alignment, payload a high','local horizonal alignment','','','')
            
            figure(2)
            NRac = obj.simout.NRac;
            NRbc = obj.simout.NRbc;
            NRdc = obj.simout.NRdc;
            NRec = obj.simout.NRec;
            NVac = obj.simout.NVac;
            NVbc = obj.simout.NVbc;
            NVdc = obj.simout.NVdc;
            NVec = obj.simout.NVec;
            for i=length(NRac):-1:1
                Hac = obj.ma*cross(NRac(i,:),NVac(i,:));
                Hbc = obj.mb*cross(NRbc(i,:),NVbc(i,:));
                Hdc = obj.md*cross(NRdc(i,:),NVdc(i,:));
                Hec = obj.me*cross(NRec(i,:),NVec(i,:));
                hc(i) = norm(Hac+Hbc+Hdc+Hec);
            end
            ax = axes();
            hold on
            plot(ax,t,hc);
            plot(ax,t(idxMinDiff),hc(idxMinDiff),'ro',LineWidth=2)
            plot(ax,t(idxMaxDiff),hc(idxMaxDiff),'bo',LineWidth=2)
            plot(ax,t(idx0Diff),hc(idx0Diff),'ko',LineWidth=2)
            title("Spacecraft angular momentum")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            legend('','local vertical alignment, a low','local vertical alignment, a high','local horizonal alignment')

            figure(3)
            axL = subplot(1,2,1);
            title("Spacecraft Attitude Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axL,3)
            data = [NRac;NRbc;NRdc;NRec];
            axis(getAxisLimits3(data))
            hold on
            plot3(axL,NRac(:,1),NRac(:,2),NRac(:,3),'k')
            plot3(axL,NRbc(:,1),NRbc(:,2),NRbc(:,3),'k')
            plot3(axL,[NRdc(end,1) NRac(end,1)],[NRdc(end,2) NRac(end,2)],[NRdc(end,3) NRac(end,3)],'b')
            plot3(axL,[NRec(end,1) NRac(end,1)],[NRec(end,2) NRac(end,2)],[NRec(end,3) NRac(end,3)],'b')
            plot3(axL,[NRdc(end,1) NRbc(end,1)],[NRdc(end,2) NRbc(end,2)],[NRdc(end,3) NRbc(end,3)],'b')
            plot3(axL,[NRec(end,1) NRbc(end,1)],[NRec(end,2) NRbc(end,2)],[NRec(end,3) NRbc(end,3)],'b')
            plot3(axL,[NRdc(end,1) NRec(end,1)],[NRdc(end,2) NRec(end,2)],[NRdc(end,3) NRec(end,3)],'k--')
            legend('Payload a & b', '', 'tethers', '', '', '', 'truss')
            axR = subplot(1,2,2);
            title("Spacecraft Orbit Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axR,3)
            NRdo = obj.simout.NRdo;
            NReo = obj.simout.NReo;
            data = [NRdo;NReo];
            axis(getAxisLimits3(data))
            hold on
            plot3(axR,NRdo(:,1),NRdo(:,2),NRdo(:,3))
            plot3(axR,NReo(:,1),NReo(:,2),NReo(:,3))
            legend('Truss end (D)','Truss end (E)')
        end

        function animate(obj)
            % t = obj.simout.tout;
            % rounded_values = floor(floor(t));
            % unique_values = unique(rounded_values);
            % idx = zeros(length(unique_values),1);
            % for i = 1:length(unique_values)
            %     index = find(rounded_values == unique_values(i), 1, 'first');
            %     idx(i) = index;
            % end
            NRac = obj.simout.NRac;
            NRbc = obj.simout.NRbc;
            NRdc = obj.simout.NRdc;
            NRec = obj.simout.NRec;
            NRco = obj.simout.NRco;
            figure
            axL = subplot(1,2,1);
            title("Spacecraft Attitude Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axL,3)
            data = [NRac;NRbc;NRdc;NRec];
            axis(getAxisLimits3(data))
            trace_a = animatedline(axL);
            trace_b = animatedline(axL);
            tether_ad = animatedline(axL);
            tether_ae = animatedline(axL);
            tether_bd = animatedline(axL);
            tether_be = animatedline(axL);
            truss_de = animatedline(axL);
            axR = subplot(1,2,2);
            title("Spacecraft Orbit Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axR,3)
            data = NRco;
            limits = getAxisLimits3(data);

            axis([min(limits) max(limits) min(limits) max(limits) min(limits) max(limits)])
            truss_center = animatedline(axR);
            pause;
            for i=1:20:length(NRac)
                addpoints(trace_a,NRac(i,1),NRac(i,2),NRac(i,3))
                addpoints(trace_b,NRbc(i,1),NRbc(i,2),NRbc(i,3))
                addpoints(tether_ad,[NRdc(i,1) NRac(i,1)],[NRdc(i,2) NRac(i,2)],[NRdc(i,3) NRac(i,3)])
                addpoints(tether_ae,[NRec(i,1) NRac(i,1)],[NRec(i,2) NRac(i,2)],[NRec(i,3) NRac(i,3)])
                addpoints(tether_bd,[NRdc(i,1) NRbc(i,1)],[NRdc(i,2) NRbc(i,2)],[NRdc(i,3) NRbc(i,3)])
                addpoints(tether_be,[NRec(i,1) NRbc(i,1)],[NRec(i,2) NRbc(i,2)],[NRec(i,3) NRbc(i,3)])
                addpoints(truss_de,[NRdc(i,1) NRec(i,1)],[NRdc(i,2) NRec(i,2)],[NRdc(i,3) NRec(i,3)])
                addpoints(truss_center,NRco(i,1),NRco(i,2),NRco(i,3))
                drawnow limitrate
                clearpoints(tether_ad)
                clearpoints(tether_ae)
                clearpoints(tether_bd)
                clearpoints(tether_be)
                clearpoints(truss_de)
            end
        end
    end
end
