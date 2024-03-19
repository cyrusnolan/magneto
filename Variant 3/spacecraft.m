classdef spacecraft
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        NX
        model
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
            obj.model = "dynamics";
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
            preloadFcn = "gravity = "+gravity+"; " + newline + ...
                "k = "+obj.p.k+"; " + newline + ...
                "b = "+obj.p.b+"; " + newline + ...
                "ktr = "+obj.p.ktr+"; " + newline + ...
                "btr = "+obj.p.btr+"; " + newline + ...
                "d0 = "+obj.p.d0+"; " + newline + ...
                "l0 = "+obj.p.l0+"; " + newline + ...
                "mu = "+obj.p.mu+"; " + newline + ...
                "mp = "+obj.p.mp+"; " + newline + ...
                "mp_inv = "+1/obj.p.mp+"; " + newline + ...
                "mtr = "+obj.p.mtr+"; " + newline + ...
                "mtr_inv = "+1/obj.p.mtr+"; " + newline + ...
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
            switch config
                case "v1"
                    set_param(obj.model+"/Mass a", "Commented", "off")
                    set_param(obj.model+"/Mass b", "Commented", "off")
                    set_param(obj.model+"/Mass f", "Commented", "on")
                    set_param(obj.model+"/Mass g", "Commented", "on")
                case "octahedron"
                    set_param(obj.model+"/Mass a", "Commented", "off")
                    set_param(obj.model+"/Mass b", "Commented", "off")
                    set_param(obj.model+"/Mass f", "Commented", "off")
                    set_param(obj.model+"/Mass g", "Commented", "off")
            end
            save_system(obj.model)
            close_system(obj.model)
            disp("Simulating configuration: " +config)
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
            % HG = angular momentum of spacecraft about its center of mass
            % HGo = angular momentum of spacecraft center of mass about Earth's center
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
                HGo(i) = totalMass*norm(cross(NRGo(i,:), NVGo(i,:)));
                HaG = mp*cross(NRaG(i,:), NVaG(i,:));
                HbG = mp*cross(NRbG(i,:), NVbG(i,:));
                HfG = mp*cross(NRfG(i,:), NVfG(i,:));
                HgG = mp*cross(NRgG(i,:), NVgG(i,:));
                HdG = mtr*cross(NRdG(i,:),NVdG(i,:));
                HeG = mtr*cross(NReG(i,:),NVeG(i,:));
                HG(i) = norm(HaG + HbG + HfG + HgG + HdG + HeG);
                Hao = mp*cross(NRao(i,:), NVao(i,:));
                Hbo = mp*cross(NRbo(i,:), NVbo(i,:));
                Hfo = mp*cross(NRfo(i,:), NVfo(i,:));
                Hgo = mp*cross(NRgo(i,:), NVgo(i,:));
                Hdo = mtr*cross(NRdo(i,:), NVdo(i,:));
                Heo = mtr*cross(NReo(i,:), NVeo(i,:));
                Ho(i) = norm(Hao + Hbo + Hfo + Hgo + Hdo + Heo);
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
            prepFigPresentation2(obj, gcf)
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
            [idx_lv_ab, idx_lh_ab] = lvlh_ab(obj);
            [idx_lv_fg, idx_lh_fg] = lvlh_fg(obj);
            [HGo, HG, Ho] = getAngularMom(obj);

            figure;
            hold on
            % plot(t, HGo);
            HGo_lv_ab = HGo(idx_lv_ab);
            t_HGo_lv_ab = t(idx_lv_ab);
            HGo_lh_ab = HGo(idx_lh_ab);
            t_HGo_lh_ab = t(idx_lh_ab);
            [idx_lv_HGo_max, ~] = localMax(HGo_lv_ab);
            [idx_lv_HGo_min, ~] = localMin(HGo_lv_ab);
            [idx_lh_HGo_max, ~] = localMax(HGo_lh_ab);
            [idx_lh_HGo_min, ~] = localMin(HGo_lh_ab);
            plot(t_HGo_lv_ab, HGo_lv_ab)
            plot(t_HGo_lh_ab, HGo_lh_ab)
            plot(t_HGo_lv_ab(idx_lv_HGo_max), HGo_lv_ab(idx_lv_HGo_max), 'k--')
            plot(t_HGo_lv_ab(idx_lv_HGo_min), HGo_lv_ab(idx_lv_HGo_min), 'k--')
            plot(t_HGo_lh_ab(idx_lh_HGo_max), HGo_lh_ab(idx_lh_HGo_max), 'k--')
            plot(t_HGo_lh_ab(idx_lh_HGo_min), HGo_lh_ab(idx_lh_HGo_min), 'k--')
            title("Orbital Angular Momentum")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            legend('ab vertical', 'ab horizonal')
            grid on;
            prepFigPresentation2(obj, gcf)

            figure;
            hold on
            %plot(t, HG)
            HG_lv_ab = HG(idx_lv_ab);
            t_HG_lv_ab = t(idx_lv_ab);
            HG_lh_ab = HG(idx_lh_ab);
            t_HG_lh_ab = t(idx_lh_ab);
            [idx_lv_HG_max, ~] = localMax(HG_lv_ab);
            [idx_lv_HG_min, ~] = localMin(HG_lv_ab);
            [idx_lh_HG_max, ~] = localMax(HG_lh_ab);
            [idx_lh_HG_min, ~] = localMin(HG_lh_ab);
            plot(t_HG_lv_ab, HG_lv_ab)
            plot(t_HG_lh_ab, HG_lh_ab)
            plot(t_HG_lv_ab(idx_lv_HG_max), HG_lv_ab(idx_lv_HG_max), 'k--')
            plot(t_HG_lv_ab(idx_lv_HG_min), HG_lv_ab(idx_lv_HG_min), 'k--')
            plot(t_HG_lh_ab(idx_lh_HG_max), HG_lh_ab(idx_lh_HG_max), 'k--')
            plot(t_HG_lh_ab(idx_lh_HG_min), HG_lh_ab(idx_lh_HG_min), 'k--')
            title("Spacecraft Angular Momentum about its Center of Mass")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            legend('ab vertical', 'ab horizonal')
            grid on;
            prepFigPresentation2(obj, gcf)

            figure;
            ax = axes;
            hold on
            plot(ax, t, Ho)
            ylim([0.95*min(Ho)  1.05*max(Ho)])
            title("Spacecraft Total Angular Momentum")
            xlabel("t (s)")
            ylabel("angular momentum (kg*m^2/s)")
            prepFigPresentation2(obj, gcf)
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
            axis(getAxisLimits3(dataR))
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

        function prepFigPresentation2(~, fignum)
        %
        % prepares a figure for presentations
        %
        % Fontsize: 14
        % Fontweight: bold
        % LineWidth: 1.5
        
            figure(fignum);
            fig_children=get(fignum,'children'); %find all sub-plots
            
            for i=1:length(fig_children)
            
                if class(fig_children(i)) == "matlab.graphics.illustration.Legend"
                    set(fig_children(i),'FontSize',16,'Interpreter','latex');
                end

                set(fig_children(i),'FontSize',16);
                set(fig_children(i),'FontWeight','bold');
                
                fig_children_children=get(fig_children(i),'Children');
                set(fig_children_children,'LineWidth',1.5);
            end
        end
    end
end
