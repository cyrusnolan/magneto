classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        NX
        simout
        model
    end
    
    methods

        function obj = spacecraft(p, NX)
            arguments
                p {mustBeA(p,"parameters")}
                NX {mustBeA(NX,"state_vec_N")}
            end
            obj.p = p;
            obj.NX = NX; 
            obj.model = p.model;
        end

        function obj = sim(obj, stopTime, gravity)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                stopTime {mustBePositive,mustBeReal,mustBeNumeric}
                gravity {mustBeA(gravity,"logical")} = true
            end
            preloadFcn = "gravity = "+gravity+"; " + newline + ...
                "k = "+obj.p.k+"; " + newline + ...
                "b = "+obj.p.b+"; " + newline + ...
                "ktr = "+obj.p.ktr+"; " + newline + ...
                "btr = "+obj.p.btr+"; " + newline + ...
                "ks = "+obj.p.ks+"; " + newline + ...
                "bs = "+obj.p.bs+"; " + newline + ...
                "d0 = "+obj.p.d0+"; " + newline + ...
                "l0 = "+obj.p.l0+"; " + newline + ...
                "ls0 = "+obj.p.ls0+"; " + newline + ...
                "mu = "+obj.p.mu+"; " + newline + ...
                "mp = "+obj.p.mp+"; " + newline + ...
                "mp_inv = "+1/obj.p.mp+"; " + newline + ...
                "mtr = "+obj.p.mtr+"; " + newline + ...
                "mtr_inv = "+1/obj.p.mtr+"; " + newline + ...
                "invIm = "+obj.p.invIm+"; " + newline + ...
                "rm = "+obj.p.rm+"; " + newline + ...
                "NRao0 = "+mat2str(obj.NX.NRao0)+"; " + newline + ...
                "NRbo0 = "+mat2str(obj.NX.NRbo0)+"; " + newline + ...
                "NRfo0 = "+mat2str(obj.NX.NRfo0)+"; " + newline + ...
                "NRgo0 = "+mat2str(obj.NX.NRgo0)+"; " + newline + ...
                "NRdo0 = "+mat2str(obj.NX.NRdo0)+"; " + newline + ...
                "NReo0 = "+mat2str(obj.NX.NReo0)+"; " + newline + ...
                "NVao0 = "+mat2str(obj.NX.NVao0)+"; " + newline + ...
                "NVbo0 = "+mat2str(obj.NX.NVbo0)+"; " + newline + ...
                "NVfo0 = "+mat2str(obj.NX.NVfo0)+"; " + newline + ...
                "NVgo0 = "+mat2str(obj.NX.NVgo0)+"; " + newline + ...
                "NVdo0 = "+mat2str(obj.NX.NVdo0)+"; " + newline + ...
                "NVeo0 = "+mat2str(obj.NX.NVeo0)+"; ";                     
            load_system(obj.model)
            set_param(obj.model, 'PreloadFcn', preloadFcn);
            postloadFcn = "StopTime = "+stopTime+"; ";
            set_param(obj.model, 'PostloadFcn', postloadFcn);
            save_system(obj.model)
            close_system(obj.model)
            disp("Simulating model: " +obj.model)
            obj.simout = sim(obj.model+".slx");
        end

        function [idx_lv, idx_lh] = lvlh_ab(obj)
            % find indicies corresponding to local vertical or local
            % horizontal configuration for mass a and b
            NRao = obj.simout.NRao;
            NRbo = obj.simout.NRbo;
            for i=length(NRao):-1:1
                rao = norm(NRao(i,:));
                rbo = norm(NRbo(i,:));
                diff(i) = rao-rbo;
            end
            diffMaxCond = [false, (diff(2:end-1) > diff(1:end-2)) & (diff(2:end-1) > diff(3:end)), false];
            diffMinCond = [false, (diff(2:end-1) < diff(1:end-2)) & (diff(2:end-1) < diff(3:end)), false];
            diff0Cond = [false, (abs(diff(2:end-1)) < abs(diff(1:end-2))) & (abs(diff(2:end-1)) < abs(diff(3:end))), false];
            idx_aLow = find(diffMinCond);
            idx_aHigh = find(diffMaxCond);
            idx_lv = sort([idx_aLow  idx_aHigh]);
            idx_lh = find(diff0Cond);
        end

        function [idx_lv, idx_lh] = lvlh_fg(obj)
            % find indicies corresponding to local vertical or local
            % horizontal configuration for mass f and g
            NRfo = obj.simout.NRfo;
            NRgo = obj.simout.NRgo;
            for i=length(NRfo):-1:1
                rfo = norm(NRfo(i,:));
                rgo = norm(NRgo(i,:));
                diff(i) = rfo-rgo;
            end
            diffMaxCond = [false, (diff(2:end-1) > diff(1:end-2)) & (diff(2:end-1) > diff(3:end)), false];
            diffMinCond = [false, (diff(2:end-1) < diff(1:end-2)) & (diff(2:end-1) < diff(3:end)), false];
            diff0Cond = [false, (abs(diff(2:end-1)) < abs(diff(1:end-2))) & (abs(diff(2:end-1)) < abs(diff(3:end))), false];
            idx_fLow = find(diffMinCond);
            idx_fHigh = find(diffMaxCond);
            idx_lv = sort([idx_fLow  idx_fHigh]);
            idx_lh = find(diff0Cond);
        end

        function [HGo, HG, Ho] = getAngularMom(obj)
            % HGo = angular momentum of spacecraft center of mass about Earth's center
            % HG = angular momentum of spacecraft about its center of mass
            % Ho = total angular momentum of spacecraft about Earth's center
            mp = obj.p.mp;
            mtr = obj.p.mtr;
            totalMass = 4*mp + 2*mtr;
            NRao = obj.simout.NRao;
            NRbo = obj.simout.NRbo;
            NRfo = obj.simout.NRfo;
            NRgo = obj.simout.NRgo;
            NRdo = obj.simout.NRdo;
            NReo = obj.simout.NReo;
            NVao = obj.simout.NVao;
            NVbo = obj.simout.NVbo;
            NVfo = obj.simout.NVfo;
            NVgo = obj.simout.NVgo;
            NVdo = obj.simout.NVdo;
            NVeo = obj.simout.NVeo;
            NRGo = (mp*(NRao + NRbo + NRfo + NRgo) + mtr*(NRdo + NReo))/totalMass;
            NVGo = (mp*(NVao + NVbo + NVfo + NVgo) + mtr*(NVdo + NVeo))/totalMass;
            NRaG = NRao - NRGo;
            NRbG = NRbo - NRGo;
            NRfG = NRfo - NRGo;
            NRgG = NRgo - NRGo;
            NRdG = NRdo - NRGo;
            NReG = NReo - NRGo;
            NVaG = NVao - NVGo;
            NVbG = NVbo - NVGo;
            NVfG = NVfo - NVGo;
            NVgG = NVgo - NVGo;
            NVdG = NVdo - NVGo;
            NVeG = NVeo - NVGo;
            for i=length(NRao):-1:1
                % HGo
                HGo(i,1:3) = totalMass*(cross(NRGo(i,:), NVGo(i,:)));
                % HG
                HaG = mp*cross(NRaG(i,:), NVaG(i,:));
                HbG = mp*cross(NRbG(i,:), NVbG(i,:));
                HfG = mp*cross(NRfG(i,:), NVfG(i,:));
                HgG = mp*cross(NRgG(i,:), NVgG(i,:));
                HdG = mtr*cross(NRdG(i,:),NVdG(i,:));
                HeG = mtr*cross(NReG(i,:),NVeG(i,:));
                HG(i,1:3) = HaG + HbG + HfG + HgG + HdG + HeG;
                % Ho
                Hao = mp*cross(NRao(i,:), NVao(i,:));
                Hbo = mp*cross(NRbo(i,:), NVbo(i,:));
                Hfo = mp*cross(NRfo(i,:), NVfo(i,:));
                Hgo = mp*cross(NRgo(i,:), NVgo(i,:));
                Hdo = mtr*cross(NRdo(i,:), NVdo(i,:));
                Heo = mtr*cross(NReo(i,:), NVeo(i,:));
                Ho(i,1:3) = Hao + Hbo + Hfo + Hgo + Hdo + Heo;
            end
        end

        function [a, e, incl, W] = getOrbitalElements(obj)
            NRao = obj.simout.NRao;
            NRbo = obj.simout.NRbo;
            NRfo = obj.simout.NRfo;
            NRgo = obj.simout.NRgo;
            NRdo = obj.simout.NRdo;
            NReo = obj.simout.NReo;
            NVao = obj.simout.NVao;
            NVbo = obj.simout.NVbo;
            NVfo = obj.simout.NVfo;
            NVgo = obj.simout.NVgo;
            NVdo = obj.simout.NVdo;
            NVeo = obj.simout.NVeo;
            mp = obj.p.mp;
            mtr = obj.p.mtr;
            totalMass = 4*mp + 2*mtr;
            Nr = (mp*(NRao + NRbo + NRfo + NRgo) + mtr*(NRdo + NReo))/totalMass;
            Nv = (mp*(NVao + NVbo + NVfo + NVgo) + mtr*(NVdo + NVeo))/totalMass;
            mu = obj.p.mu;
            for i = length(Nr):-1:1
                % angular momentum
                Nh = cross(Nr(i,:), Nv(i,:));
                h = norm(Nh);
                % inclination
                incl(i) = acos(Nh(3)/h);
                % eccentricity
                rhat = Nr(i,:)/norm(Nr(i,:));
                Ne = cross(Nv(i,:), Nh)/mu - rhat;
                e(i) = norm(Ne);
                % RAAN
                Nn = cross([0 0 1], Nh);
                n = norm(Nn);
                W(i) = acos(Nn(1)/n);
                % semimajor axis
                a(i) = h^2/mu/(1-e(i)^2);
            end
        end

        function plotTrussLength(obj)
            t = obj.simout.tout;
            rde = obj.simout.rde;
            [idx_lv_ab, idx_lh_ab] = lvlh_ab(obj);
            [idx_lv_fg, idx_lh_fg] = lvlh_fg(obj);

            figure;
            hold on
            plot(t, rde)
            plot(t(idx_lv_ab), rde(idx_lv_ab), 'ro')
            plot(t(idx_lh_ab), rde(idx_lh_ab), 'ko')
            plot(t(idx_lv_fg), rde(idx_lv_fg), 'bo')
            plot(t(idx_lh_fg), rde(idx_lh_fg), 'mo')
            title("Truss Length vs time")
            xlabel("t (s)")
            ylabel("truss length (m)")
            legend('', 'ab vertical', 'ab horizonal', 'fg vertical', 'fg horizontal')
            prepFigPresentation2(gcf)
        end

        function plotTetherLength(obj)
            t = obj.simout.tout;
            [idx_lv_ab, idx_lh_ab] = lvlh_ab(obj);
            [idx_lv_fg, idx_lh_fg] = lvlh_fg(obj);
            figure;

            ax11 = subplot(2, 2, 1);
            hold on
            plot(ax11, t, obj.simout.rad, 'k')
            plot(ax11, t, obj.simout.rae, 'r-.')
            plot(ax11, t(idx_lv_ab), obj.simout.rad(idx_lv_ab), 'ro')
            plot(ax11, t(idx_lh_ab), obj.simout.rad(idx_lh_ab), 'ko')
            title("mass a")
            xlabel("t (s)")
            ylabel("tether length (m)")
            legend("tether one","tether two","local vertical","local horizontal")
            
            ax12 = subplot(2,2,2);
            hold on
            plot(ax12, t, obj.simout.rbd, 'k')
            plot(ax12, t, obj.simout.rbe, 'r-.')
            plot(ax12, t(idx_lv_ab), obj.simout.rbd(idx_lv_ab), 'ro')
            plot(ax12, t(idx_lh_ab), obj.simout.rbd(idx_lh_ab), 'ko')
            title("mass b")
            xlabel("t (s)")
            ylabel("tether length (m)")
            legend("tether one","tether two","local vertical","local horizontal")

            ax21 = subplot(2,2,3);
            hold on
            plot(ax21, t, obj.simout.rfd, 'k')
            plot(ax21, t, obj.simout.rfe, 'r-.')
            plot(ax21, t(idx_lv_fg), obj.simout.rfd(idx_lv_fg), 'ro')
            plot(ax21, t(idx_lh_fg), obj.simout.rfd(idx_lh_fg), 'ko')
            title("mass f")
            xlabel("t (s)")
            ylabel("tether length (m)")
            legend("tether one","tether two","local vertical","local horizontal")

            ax22 = subplot(2,2,4);
            hold on
            plot(ax22, t, obj.simout.rgd, 'k')
            plot(ax22, t, obj.simout.rge, 'r-.')
            plot(ax22, t(idx_lv_fg), obj.simout.rgd(idx_lv_fg), 'ro')
            plot(ax22, t(idx_lh_fg), obj.simout.rgd(idx_lh_fg), 'ko')
            title("mass g")
            xlabel("t (s)")
            ylabel("tether length (m)")
            legend("tether one","tether two","local vertical","local horizontal")
        end

        function plotAngularMom(obj)
            t = obj.simout.tout;
            [HGo, HG, Ho] = getAngularMom(obj);
            
            figure;
            sgtitle("Orbital Angular Momentum")
            for i = 3:-1:1
                spf1(i) = subplot(3,1,i);
                plot(t, HGo(:,i))
                grid on
                xlabel("time (s)")
                ylabel("kg*m^2/s");
                if i == 1
                    title("ECI X");
                elseif i == 2
                    title("ECI Y");
                elseif i == 3
                    title("ECI Z");
                end
            end
            prepFigPresentation2(gcf)
            linkaxes_y(spf1);

            figure;
            sgtitle("Spacecraft Spin Angular Momentum")
            for i = 3:-1:1
                spf2(i) = subplot(3,1,i);
                plot(t, HG(:,i))
                grid on
                xlabel("time (s)")
                ylabel("kg*m^2/s");
                if i == 1
                    title("ECI X");
                elseif i == 2
                    title("ECI Y");
                elseif i == 3
                    title("ECI Z");
                end
            end
            prepFigPresentation2(gcf)
            linkaxes_y(spf2);

            figure;
            sgtitle("Total Angular Momentum")
            for i = 3:-1:1
                spf3(i) = subplot(3,1,i);
                plot(t, Ho(:,i))
                grid on
                xlabel("time (s)")
                ylabel("kg*m^2/s");
                %spf3(i).YAxis.TickLabelFormat = '%.1f';
                if i == 1
                    title("ECI X");
                elseif i == 2
                    title("ECI Y");
                elseif i == 3
                    title("ECI Z");
                end
            end
            prepFigPresentation2(gcf)
            linkaxes_y(spf3);
        end

        function plotOrbitalElements(obj)
            [a, e, incl, W] = getOrbitalElements(obj);
            t = obj.simout.tout;

            figure;
            plot(t, a)
            title("Semi-Major Axis")
            xlabel("t (s)")
            ylabel("(m)")
            prepFigPresentation2(gcf)

            figure;
            plot(t, e)
            title("Eccentricity")
            xlabel("t (s)")
            prepFigPresentation2(gcf)

            figure;
            subplot(1,2,1)
            plot(t, rad2deg(incl))
            title("Inclination")
            xlabel("t (s)")
            ylabel("deg")

            subplot(1,2,2)
            plot(t, rad2deg(W))
            title("Longitude of Ascending Node")
            xlabel("t (s)")
            ylabel("deg")
            prepFigPresentation2(gcf)
        end

        function animate(obj)
            NRac = obj.simout.NRac;
            NRbc = obj.simout.NRbc;
            NRfc = obj.simout.NRfc;
            NRgc = obj.simout.NRgc;
            NRdc = obj.simout.NRdc;
            NRec = obj.simout.NRec;
            NRco = obj.simout.NRco;

            figure;
            axL = subplot(1,2,1);
            axis equal
            title("Spacecraft Attitude Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axL,3)
            dataL = [NRac;NRbc;NRfc;NRgc;NRdc;NRec];
            axis(getAxisLimits3(dataL))
            trace_a = animatedline(axL);
            trace_b = animatedline(axL);
            tether_ad = animatedline(axL);
            tether_ae = animatedline(axL);
            tether_bd = animatedline(axL);
            tether_be = animatedline(axL);
            trace_f = animatedline(axL);
            trace_g = animatedline(axL);
            tether_fd = animatedline(axL);
            tether_fe = animatedline(axL);
            tether_gd = animatedline(axL);
            tether_ge = animatedline(axL);
            truss_de = animatedline(axL);

            axR = subplot(1,2,2);
            title("Spacecraft Orbit Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axR,3)
            dataR = NRco;
            lims = getAxisLimits3(dataR);
            axis([lims(1:2) lims(1:2) lims(5:6)]);
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
                addpoints(trace_f,NRfc(i,1),NRfc(i,2),NRfc(i,3))
                addpoints(trace_g,NRgc(i,1),NRgc(i,2),NRgc(i,3))
                addpoints(tether_fd,[NRdc(i,1) NRfc(i,1)],[NRdc(i,2) NRfc(i,2)],[NRdc(i,3) NRfc(i,3)])
                addpoints(tether_fe,[NRec(i,1) NRfc(i,1)],[NRec(i,2) NRfc(i,2)],[NRec(i,3) NRfc(i,3)])
                addpoints(tether_gd,[NRdc(i,1) NRgc(i,1)],[NRdc(i,2) NRgc(i,2)],[NRdc(i,3) NRgc(i,3)])
                addpoints(tether_ge,[NRec(i,1) NRgc(i,1)],[NRec(i,2) NRgc(i,2)],[NRec(i,3) NRgc(i,3)])
                addpoints(truss_center,NRco(i,1),NRco(i,2),NRco(i,3))
                drawnow limitrate
                clearpoints(tether_ad)
                clearpoints(tether_ae)
                clearpoints(tether_bd)
                clearpoints(tether_be)
                clearpoints(truss_de)
                clearpoints(tether_fd)
                clearpoints(tether_fe)
                clearpoints(tether_gd)
                clearpoints(tether_ge)
            end
            addpoints(tether_ad,[NRdc(end,1) NRac(end,1)],[NRdc(end,2) NRac(end,2)],[NRdc(end,3) NRac(end,3)])
            addpoints(tether_ae,[NRec(end,1) NRac(end,1)],[NRec(end,2) NRac(end,2)],[NRec(end,3) NRac(end,3)])
            addpoints(tether_bd,[NRdc(end,1) NRbc(end,1)],[NRdc(end,2) NRbc(end,2)],[NRdc(end,3) NRbc(end,3)])
            addpoints(tether_be,[NRec(end,1) NRbc(end,1)],[NRec(end,2) NRbc(end,2)],[NRec(end,3) NRbc(end,3)])
            addpoints(truss_de,[NRdc(end,1) NRec(end,1)],[NRdc(end,2) NRec(end,2)],[NRdc(end,3) NRec(end,3)])
            addpoints(tether_fd,[NRdc(end,1) NRfc(end,1)],[NRdc(end,2) NRfc(end,2)],[NRdc(end,3) NRfc(end,3)])
            addpoints(tether_fe,[NRec(end,1) NRfc(end,1)],[NRec(end,2) NRfc(end,2)],[NRec(end,3) NRfc(end,3)])
            addpoints(tether_gd,[NRdc(end,1) NRgc(end,1)],[NRdc(end,2) NRgc(end,2)],[NRdc(end,3) NRgc(end,3)])
            addpoints(tether_ge,[NRec(end,1) NRgc(end,1)],[NRec(end,2) NRgc(end,2)],[NRec(end,3) NRgc(end,3)])
        end

    end
end
