classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        NX
        simout
    end
    
    methods
        function obj = spacecraft(p, NX)
            arguments
                p {mustBeA(p,"parameters")}
                NX {mustBeA(NX,"state_vec_N")}
            end
            obj.p = p;
            obj.NX = NX; 
        end

        function obj = sim(obj, stopTime, gravity)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                stopTime {mustBePositive,mustBeReal,mustBeNumeric}
                gravity {mustBeA(gravity,"logical")} = true
            end
            config = obj.p.config;
            model = "dynamics";
            varInit = "StopTime = "+stopTime+"; " + ...
                "gravity = "+gravity+"; " + ...
                "k = "+obj.p.k+"; " + ...
                "b = "+obj.p.b+"; " + ...
                "ktr = "+obj.p.ktr+"; " + ...
                "btr = "+obj.p.btr+"; " + ...
                "d0 = "+obj.p.d0+"; " + ...
                "l0 = "+obj.p.l0+"; " + ...
                "mu = "+obj.p.mu+"; " + ...
                "mp = "+obj.p.mp+"; " + ...
                "mp_inv = "+1/obj.p.mp+"; " + ...
                "mtr = "+obj.p.mtr+"; " + ...
                "mtr_inv = "+1/obj.p.mtr+"; " + ...
                "NRao0 = "+mat2str(obj.NX.NRao0)+"; " + ...
                "NRbo0 = "+mat2str(obj.NX.NRbo0)+"; " + ...
                "NRdo0 = "+mat2str(obj.NX.NRdo0)+"; " + ...
                "NReo0 = "+mat2str(obj.NX.NReo0)+"; " + ...
                "NVao0 = "+mat2str(obj.NX.NVao0)+"; " + ...
                "NVbo0 = "+mat2str(obj.NX.NVbo0)+"; " + ...
                "NVdo0 = "+mat2str(obj.NX.NVdo0)+"; " + ...
                "NVeo0 = "+mat2str(obj.NX.NVeo0)+"; ";                     
            load_system(model)
            switch config
                case "v1"
                    set_param(model+"/Mass a", "Commented", "off")
                    set_param(model+"/Mass b", "Commented", "off")
                    set_param(model+"/Mass f", "Commented", "on")
                    set_param(model+"/Mass g", "Commented", "on")
                case "octahedron"
                    set_param(model+"/Mass a", "Commented", "off")
                    set_param(model+"/Mass b", "Commented", "off")
                    set_param(model+"/Mass f", "Commented", "off")
                    set_param(model+"/Mass g", "Commented", "off")
                    varInit = varInit + ...
                    "NRfo0 = "+mat2str(obj.NX.NRfo0)+"; " + ...
                    "NRgo0 = "+mat2str(obj.NX.NRgo0)+"; " + ...
                    "NVfo0 = "+mat2str(obj.NX.NVfo0)+"; " + ...
                    "NVgo0 = "+mat2str(obj.NX.NVgo0)+"; ";
            end
            set_param(model, 'PreloadFcn', varInit)
            save_system(model)
            close_system(model)
            disp("Simulating configuration: " +config)
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

            if obj.p.config == "octahedron"
                NRfo = obj.simout.NRfo;
                NRgo = obj.simout.NRgo;
                for i=length(NRfo):-1:1
                    rfo = norm(NRfo(i,:));
                    rgo = norm(NRgo(i,:));
                    diff2(i) = rfo-rgo;
                end
                diffMinCond2 = [false, (diff2(2:end-1) > diff2(1:end-2)) & (diff2(2:end-1) > diff2(3:end)), false];
                diffMaxCond2 = [false, (diff2(2:end-1) < diff2(1:end-2)) & (diff2(2:end-1) < diff2(3:end)), false];
                diff0Cond2 = [false, (abs(diff2(2:end-1)) < abs(diff2(1:end-2))) & (abs(diff2(2:end-1)) < abs(diff2(3:end))), false];
                idxMinDiff2 = find(diffMinCond2);
                idxMaxDiff2 = find(diffMaxCond2);
                idx0Diff2 = find(diff0Cond2);
            end

            figure() % truss length (axes=axL) and tether length (axes=axR) vs time
            axL = subplot(1,2,1);
            hold on
            t = obj.simout.tout;
            rde = obj.simout.rde;
            plot(axL,t,rde)
            plot(axL,t(idxMinDiff),rde(idxMinDiff),'ro',LineWidth=1)
            plot(axL,t(idxMaxDiff),rde(idxMaxDiff),'ro',LineWidth=1)
            plot(axL,t(idx0Diff),rde(idx0Diff),'ko',LineWidth=1)
            if obj.p.config == "octahedron"
                plot(axL,t(idxMinDiff2),rde(idxMinDiff2),'bo',LineWidth=1)
                plot(axL,t(idxMaxDiff2),rde(idxMaxDiff2),'bo',LineWidth=1)
                plot(axL,t(idx0Diff2),rde(idx0Diff2),'mo',LineWidth=1)
            end
            title("Distance between mass D and E (Truss Length) vs time")
            xlabel("t (s)")
            ylabel("truss length (m)")
            if obj.p.config == "v1"
                legend('','local vertical','','local horizonal')
            elseif obj.p.config == "octahedron"
                legend('','ab vertical','','ab horizonal','fg vertical','','fg horizontal')
            end

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
                if obj.p.config == "octahedron"
                    NRfc = obj.simout.NRfc;
                    NRgc = obj.simout.NRgc;
                    NVfc = obj.simout.NVfc;
                    NVgc = obj.simout.NVgc;
                    Hfc = obj.p.mp*cross(NRfc(i,:),NVfc(i,:));
                    Hgc = obj.p.mp*cross(NRgc(i,:),NVgc(i,:));
                    hc(i) = norm(Hac+Hbc+Hdc+Hec+Hfc+Hgc);
                else
                    hc(i) = norm(Hac+Hbc+Hdc+Hec);
                end
            end
            axL = subplot(1,2,1);
            hold on
            plot(axL,t,hc);
            plot(axL,t(idxMinDiff),hc(idxMinDiff),'ro',LineWidth=2)
            plot(axL,t(idxMaxDiff),hc(idxMaxDiff),'ro',LineWidth=2)
            plot(axL,t(idx0Diff),hc(idx0Diff),'ko',LineWidth=2)
            if obj.p.config == "octahedron"
                plot(axL,t(idxMinDiff2),hc(idxMinDiff2),'bo',LineWidth=1)
                plot(axL,t(idxMaxDiff2),hc(idxMaxDiff2),'bo',LineWidth=1)
                plot(axL,t(idx0Diff2),hc(idx0Diff2),'mo',LineWidth=1)
            end
            title("Spacecraft angular momentum")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            if obj.p.config == "v1"
                legend('','local vertical','','local horizonal')
            elseif obj.p.config == "octahedron"
                legend('','ab vertical','','ab horizonal','fg vertical','','fg horizontal')
            end
            axR = subplot(1,2,2);
            plot(axR,t(idx0Diff),hc(idx0Diff),'k',LineWidth=2)
            title("Spacecraft angular momentum at each ab horizontal")
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
            if obj.p.config == "octahedron"
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
            if obj.p.config == "v1"
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
            if obj.p.config == "octahedron"
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
                if obj.p.config == "octahedron"
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
                if obj.p.config == "octahedron"
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
            if obj.p.config == "octahedron"
                addpoints(tether_fd,[NRdc(end,1) NRfc(end,1)],[NRdc(end,2) NRfc(end,2)],[NRdc(end,3) NRfc(end,3)])
                addpoints(tether_fe,[NRec(end,1) NRfc(end,1)],[NRec(end,2) NRfc(end,2)],[NRec(end,3) NRfc(end,3)])
                addpoints(tether_gd,[NRdc(end,1) NRgc(end,1)],[NRdc(end,2) NRgc(end,2)],[NRdc(end,3) NRgc(end,3)])
                addpoints(tether_ge,[NRec(end,1) NRgc(end,1)],[NRec(end,2) NRgc(end,2)],[NRec(end,3) NRgc(end,3)])
            end
        end
    end
end
