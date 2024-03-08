classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p                           % parameters
        BX                          % state vector in B coordinates
        simout
    end
    
    methods
        function obj = spacecraft(parameters, BX)
            obj.p = parameters;
            obj.BX = BX;
        end

        function obj = sim(obj, stopTime, orbit, gravity)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                stopTime {mustBePositive,mustBeReal,mustBeNumeric}
                orbit {mustBeA(orbit,"orbit")}
                gravity {mustBeA(gravity,"logical")} = true
            end
            model = obj.p.model;
            NRco0 = orbit.NR;
            NVco0 = orbit.NV;
            Nh_orbit = orbit.Nh;
            Nb2 = -NRco0/norm(NRco0);
            Nb3 = Nh_orbit/norm(Nh_orbit);
            Nb1 = cross(Nb2,Nb3);
            NQB = [Nb1 Nb2 Nb3];
            NRac0 = NQB*obj.BX.BRac0;
            NRbc0 = NQB*obj.BX.BRbc0;
            NRfc0 = NQB*obj.BX.BRfc0;
            NRgc0 = NQB*obj.BX.BRgc0;
            NRdc0 = NQB*obj.BX.BRdc0;
            NRec0 = NQB*obj.BX.BRec0;
            NVac0 = NQB*obj.BX.BVac0;
            NVbc0 = NQB*obj.BX.BVbc0;
            NVfc0 = NQB*obj.BX.BVfc0;
            NVgc0 = NQB*obj.BX.BVgc0;
            NVdc0 = NQB*obj.BX.BVdc0;
            NVec0 = NQB*obj.BX.BVec0;
            NRao0 = NRac0+NRco0;
            NRbo0 = NRbc0+NRco0;
            NRfo0 = NRfc0+NRco0;
            NRgo0 = NRgc0+NRco0;
            NRdo0 = NRdc0+NRco0;
            NReo0 = NRec0+NRco0;
            NVao0 = NVac0+NVco0;
            NVbo0 = NVbc0+NVco0;
            NVfo0 = NVfc0+NVco0;
            NVgo0 = NVgc0+NVco0;
            NVdo0 = NVdc0+NVco0;
            NVeo0 = NVec0+NVco0;
            varInit = "StopTime = "+stopTime+"; " + ...
                "gravity = "+gravity+"; " + ...
                "k = "+obj.p.k+"; " + ...
                "b = "+obj.p.b+"; " + ...
                "kc = "+obj.p.kc+"; " + ...
                "bc = "+obj.p.bc+"; " + ...
                "d0 = "+obj.p.d0+"; " + ...
                "l0 = "+obj.p.l0+"; " + ...
                "mu = "+obj.p.mu+"; " + ...
                "mp = "+obj.p.mp+"; " + ...
                "mp_inv = "+1/obj.p.mp+"; " + ...
                "mtr = "+obj.p.mtr+"; " + ...
                "mtr_inv = "+1/obj.p.mtr+"; " + ...
                "NRao0 = "+mat2str(NRao0)+"; " + ...
                "NRbo0 = "+mat2str(NRbo0)+"; " + ...
                "NRfo0 = "+mat2str(NRfo0)+"; " + ...
                "NRgo0 = "+mat2str(NRgo0)+"; " + ...
                "NRdo0 = "+mat2str(NRdo0)+"; " + ...
                "NReo0 = "+mat2str(NReo0)+"; " + ...
                "NVao0 = "+mat2str(NVao0)+"; " + ...
                "NVbo0 = "+mat2str(NVbo0)+"; " + ...
                "NVfo0 = "+mat2str(NVfo0)+"; " + ...
                "NVgo0 = "+mat2str(NVgo0)+"; " + ...
                "NVdo0 = "+mat2str(NVdo0)+"; " + ...
                "NVeo0 = "+mat2str(NVeo0)+";";                     
            load_system(model)
            set_param(model, 'PreloadFcn', varInit)
            save_system(model)
            close_system(model)
            disp("Simulating model: " +model)
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
                Hac = obj.p.mp*cross(NRac(i,:),NVac(i,:));
                Hbc = obj.p.mp*cross(NRbc(i,:),NVbc(i,:));
                Hdc = obj.p.mtr*cross(NRdc(i,:),NVdc(i,:));
                Hec = obj.p.mtr*cross(NRec(i,:),NVec(i,:));
                hc(i) = norm(Hac+Hbc+Hdc+Hec);
            end
            axL = subplot(1,2,1);
            hold on
            plot(axL,t,hc);
            plot(axL,t(idxMinDiff),hc(idxMinDiff),'ro',LineWidth=2)
            plot(axL,t(idxMaxDiff),hc(idxMaxDiff),'bo',LineWidth=2)
            plot(axL,t(idx0Diff),hc(idx0Diff),'ko',LineWidth=2)
            title("Spacecraft angular momentum")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            legend('','local vertical alignment, a low','local vertical alignment, a high','local horizonal alignment')
            axR = subplot(1,2,2);
            plot(axR,t(idx0Diff),hc(idx0Diff),'k',LineWidth=2)
            title("Spacecraft angular momentum at local horizontal")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")

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
            NRac = obj.simout.NRac;
            NRbc = obj.simout.NRbc;
            NRdc = obj.simout.NRdc;
            NRec = obj.simout.NRec;
            NRco = obj.simout.NRco;
            if obj.p.model == "dynamics_v2"
                NRfc = obj.simout.NRfc;
                NRgc = obj.simout.NRgc;
            end
            figure
            axL = subplot(1,2,1);
            axis equal
            title("Spacecraft Attitude Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axL,3)
            if obj.p.model == "dynamics"
                data = [NRac;NRbc;NRdc;NRec];
                axis(getAxisLimits3(data))
            else
                data = [NRac;NRbc;NRfc;NRgc;NRdc;NRec];
                axis(getAxisLimits3(data))
            end
            trace_a = animatedline(axL);
            trace_b = animatedline(axL);
            tether_ad = animatedline(axL);
            tether_ae = animatedline(axL);
            tether_bd = animatedline(axL);
            tether_be = animatedline(axL);
            truss_de = animatedline(axL);
            if obj.p.model == "dynamics_v2"
                trace_f = animatedline(axL);
                trace_g = animatedline(axL);
                tether_fd = animatedline(axL);
                tether_fe = animatedline(axL);
                tether_gd = animatedline(axL);
                tether_ge = animatedline(axL);
            end
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
                if obj.p.model == "dynamics_v2"
                    addpoints(trace_f,NRfc(i,1),NRfc(i,2),NRfc(i,3))
                    addpoints(trace_g,NRgc(i,1),NRgc(i,2),NRgc(i,3))
                    addpoints(tether_fd,[NRdc(i,1) NRfc(i,1)],[NRdc(i,2) NRfc(i,2)],[NRdc(i,3) NRfc(i,3)])
                    addpoints(tether_fe,[NRec(i,1) NRfc(i,1)],[NRec(i,2) NRfc(i,2)],[NRec(i,3) NRfc(i,3)])
                    addpoints(tether_gd,[NRdc(i,1) NRgc(i,1)],[NRdc(i,2) NRgc(i,2)],[NRdc(i,3) NRgc(i,3)])
                    addpoints(tether_ge,[NRec(i,1) NRgc(i,1)],[NRec(i,2) NRgc(i,2)],[NRec(i,3) NRgc(i,3)])
                end
                drawnow limitrate
                clearpoints(tether_ad)
                clearpoints(tether_ae)
                clearpoints(tether_bd)
                clearpoints(tether_be)
                clearpoints(truss_de)
                if obj.p.model == "dynamics_v2"
                    clearpoints(tether_fd)
                    clearpoints(tether_fe)
                    clearpoints(tether_gd)
                    clearpoints(tether_ge)
                end
            end
            addpoints(tether_ad,[NRdc(end,1) NRac(end,1)],[NRdc(end,2) NRac(end,2)],[NRdc(end,3) NRac(end,3)])
            addpoints(tether_ae,[NRec(end,1) NRac(end,1)],[NRec(end,2) NRac(end,2)],[NRec(end,3) NRac(end,3)])
            addpoints(tether_bd,[NRdc(end,1) NRbc(end,1)],[NRdc(end,2) NRbc(end,2)],[NRdc(end,3) NRbc(end,3)])
            addpoints(tether_be,[NRec(end,1) NRbc(end,1)],[NRec(end,2) NRbc(end,2)],[NRec(end,3) NRbc(end,3)])
            addpoints(truss_de,[NRdc(end,1) NRec(end,1)],[NRdc(end,2) NRec(end,2)],[NRdc(end,3) NRec(end,3)])
            if obj.p.model == "dynamics_v2"
                addpoints(tether_fd,[NRdc(end,1) NRfc(end,1)],[NRdc(end,2) NRfc(end,2)],[NRdc(end,3) NRfc(end,3)])
                addpoints(tether_fe,[NRec(end,1) NRfc(end,1)],[NRec(end,2) NRfc(end,2)],[NRec(end,3) NRfc(end,3)])
                addpoints(tether_gd,[NRdc(end,1) NRgc(end,1)],[NRdc(end,2) NRgc(end,2)],[NRdc(end,3) NRgc(end,3)])
                addpoints(tether_ge,[NRec(end,1) NRgc(end,1)],[NRec(end,2) NRgc(end,2)],[NRec(end,3) NRgc(end,3)])
            end
        end
    end
end
