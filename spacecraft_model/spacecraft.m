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

        function obj = sim(obj, model, stopTime, gravity)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            arguments
                obj
                model
                stopTime {mustBePositive,mustBeReal,mustBeNumeric}
                gravity {mustBeA(gravity,"logical")} = true
            end
            preloadFcn = "gravity = "+gravity+"; " + newline + ...
                "e_a = "+obj.p.ellipse_a+"; " + newline + ...
                "e_b = "+obj.p.ellipse_b  +"; " + newline + ...
                "k = "+obj.p.k+"; " + newline + ...
                "b = "+obj.p.b+"; " + newline + ...
                "ktr = "+obj.p.ktr+"; " + newline + ...
                "btr = "+obj.p.btr+"; " + newline + ...
                "ks = "+obj.p.ks+"; " + newline + ...
                "bs = "+obj.p.bs+"; " + newline + ...
                "d0 = "+obj.p.d0+"; " + newline + ...
                "l0 = "+obj.p.l0+"; " + newline + ...
                "p0 = "+obj.p.p0+"; " + newline + ...
                "ls0 = "+obj.p.ls0+"; " + newline + ...
                "w = "+obj.p.w+"; " + newline + ...
                "kp_amc = "+obj.p.kp_amc+"; " + newline + ...
                "ki_amc = "+obj.p.ki_amc+"; " + newline + ...
                "amc_on_time = "+obj.p.amc_on_time+"; " + newline + ...
                "kp_tl = "+obj.p.kp_tl+"; " + newline + ...
                "kd_tl = "+obj.p.kd_tl+"; " + newline + ...
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
            load_system(model)
            set_param(model, 'PreloadFcn', preloadFcn);
            postloadFcn = "StopTime = "+stopTime+"; ";
            set_param(model, 'PostloadFcn', postloadFcn);
            save_system(model)
            close_system(model)
            disp("Simulating model: " +model)
            obj.simout = sim(model+".slx");
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

        function [NRGo, NVGo] = getCG(obj)
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
            NRGo = (mp*(NRao + NRbo + NRfo + NRgo) + mtr*(NRdo + NReo))/totalMass;
            NVGo = (mp*(NVao + NVbo + NVfo + NVgo) + mtr*(NVdo + NVeo))/totalMass;
        end

        function [eap, ebp] = getTrackingError(obj)
            ea = obj.simout.error_a;
            la = obj.simout.rad;
            eap = 100*ea./la;
            eb = obj.simout.error_b;
            lb = obj.simout.rbd;
            ebp = 100*eb./lb;

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
            [NRGo, NVGo] = getCG(obj);
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
            [Nr, Nv] = getCG(obj);
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

            figure;
            hold on
            plot(t, rde)
            title("Truss Length vs time")
            xlabel("t (s)")
            ylabel("truss length (m)")
            prepFigPresentation2(gcf)
        end

        function plotTetherLength(obj, singleFig)
            arguments
                obj
                singleFig logical=true
            end
            t = obj.simout.tout;
            
            if singleFig
                figure;
                ax11 = subplot(2, 2, 1);
                hold on
                grid on
                ax12 = subplot(2,2,2);
                hold on 
                grid on
                ax21 = subplot(2,2,3);
                hold on
                grid on
                ax22 = subplot(2,2,4);
                hold on
                grid on
            else
                figure(1);
                ax11 = axes;
                hold on
                grid on
                figure(2);
                ax12 = axes;
                hold on
                figure(3);
                ax21 = axes;
                hold on
                figure(4);
                ax22 = axes;
                hold on
            end

            % mass a
            plot(ax11, t, obj.simout.rad, 'k')
            plot(ax11, t, obj.simout.rae, 'r--')
            title(ax11, "module a")
            xlabel(ax11, "time (hours)")
            ylabel(ax11, "tether length (m)")
            legend(ax11, "tether one","tether two")
            ax11.XTick = 0:1/4*60*60:max(ax11.XLim);
            ax11.XTickLabel = ax11.XTick/(60*60);
            xlim(ax11,[min(t), max(t)]);
            
            % mass b
            plot(ax12, t, obj.simout.rbd, 'k')
            plot(ax12, t, obj.simout.rbe, 'r--')
            title(ax12, "module b")
            xlabel(ax12, "time (hours)")
            ylabel(ax12, "tether length (m)")
            legend(ax12, "tether one","tether two")
            ax12.XTick = 0:1/4*60*60:max(ax12.XLim);
            ax12.XTickLabel = ax12.XTick/(60*60);
            xlim(ax12,[min(t), max(t)]);

            % mass m
            plot(ax21, t, obj.simout.rfd, 'k')
            plot(ax21, t, obj.simout.rfe, 'r--')
            title(ax21, "module m")
            xlabel(ax21, "time (hours)")
            ylabel(ax21, "tether length (m)")
            legend(ax21, "tether one","tether two")
            ax21.XTick = 0:1/4*60*60:max(ax21.XLim);
            ax21.XTickLabel = ax21.XTick/(60*60);
            xlim(ax21,[min(t), max(t)]);

            % mass n
            plot(ax22, t, obj.simout.rgd, 'k')
            plot(ax22, t, obj.simout.rge, 'r--')
            title(ax22, "module n")
            xlabel(ax22, "time (hours)")
            ylabel(ax22, "Tether Length (m)")
            legend(ax22, "tether one","tether two")
            ax22.XTick = 0:1/4*60*60:max(ax22.XLim);
            ax22.XTickLabel = ax22.XTick/(60*60);
            xlim(ax22,[min(t), max(t)]);

            prepFigPresentation(gcf)
        end

        function plotTrackingError(obj)
            [eap, ebp] = getTrackingError(obj);
            t = obj.simout.tout;

            figure;
            subplot(2,1,1)
            plot(t,eap);
            title("Module A Tracking Error")
            ylabel("percent error")
            ax = gca;
            ax.XTick = 0:1/4*60*60:max(ax.XLim);
            ax.XTickLabel = ax.XTick/(60*60);

            subplot(2,1,2)
            plot(t,ebp);
            title("Module B Tracking Error")
            xlabel("time (hours)")
            ylabel("percent error")
            ax = gca;
            ax.XTick = 0:1/4*60*60:max(ax.XLim);
            ax.XTickLabel = ax.XTick/(60*60);
            prepFigPresentation(gcf);

        end

        function plotAngularMom(obj)
            t = obj.simout.tout;
            [HGo, HG, Ho] = getAngularMom(obj);
            total_mass = 4*obj.p.mp + 2*obj.p.mtr;

            figure;
            plot(t, HG(:,2))
            title("Spin Angular Momentum")
            grid on
            xlabel("time (hours)")
            ylabel("$kg*m^2/s$");
            ax = gca;
            ax.XTick = 0:1/4*60*60:max(ax.XLim);
            ax.XTickLabel = ax.XTick/(60*60);
            xlim([min(t), max(t)]);
            linkaxes_y(ax)
            prepFigPresentation(gcf)
            ax.InnerPosition = [0.17,0.11,0.75,0.815];
            
            
            figure;
            plot(t, HGo(:,2)./total_mass)
            title("Orbital Angular Momentum")
            grid on
            xlabel("time (hours)")
            ylabel("$m^2/s$");
            ax = gca;
            ax.XTick = 0:1/4*60*60:max(ax.XLim);
            ax.XTickLabel = ax.XTick/(60*60);
            xlim([min(t), max(t)]);
            linkaxes_y(ax)
            prepFigPresentation(gcf)
            ax.InnerPosition = [0.17,0.11,0.75,0.815];

        %     figure;
        %     sgtitle("Change in Orbital Angular Momentum")
        %     for i = 3:-1:1
        %         spf1(i) = subplot(3,1,i);
        %         plot(t, HGo(:,i))
        %         grid on
        %         xlabel("time (s)")
        %         ylabel("kg*m^2/s");
        %         if i == 1
        %             text(0.9, 0.1, "ECI X", 'Units', 'normalized', 'FontSize', 16)
        %             spf1(i).XTickLabels = {};
        %             xlim([min(t), max(t)]);
        %         elseif i == 2
        %             text(0.9, 0.1, "ECI Y", 'Units', 'normalized', 'FontSize', 16)
        %             spf1(i).XTickLabels = {};
        %             xlim([min(t), max(t)]);
        %         elseif i == 3
        %             text(0.9, 0.1, "ECI Z", 'Units', 'normalized', 'FontSize', 16)
        %             xlabel(spf1(i), "time (hours)");
        %             spf1(i).XTick = 0:1/4*60*60:max(spf1(i).XLim);
        %             spf1(i).XTickLabel = spf1(i).XTick/(60*60);
        %             xlim([min(t), max(t)]);
        %         end
        %     end
        %     linkaxes_y(spf1);
        %     prepFigPresentation2(gcf)
        % 
        %     figure;
        %     sgtitle("Spacecraft Spin Angular Momentum")
        %     for i = 3:-1:1
        %         spf2(i) = subplot(3,1,i);
        %         plot(t, HG(:,i))
        %         grid on
        %         ylabel("kg*m^2/s");
        %         if i == 1
        %             text(0.9, 0.1, "ECI X", 'Units', 'normalized', 'FontSize', 16)
        %             spf2(i).XTickLabels = {};
        %             xlim([min(t), max(t)]);
        %         elseif i == 2
        %             text(0.9, 0.1, "ECI Y", 'Units', 'normalized', 'FontSize', 16)
        %             spf2(i).XTickLabels = {};
        %             xlim([min(t), max(t)]);
        %         elseif i == 3
        %             text(0.9, 0.1, "ECI Z", 'Units', 'normalized', 'FontSize', 16)
        %             xlabel(spf2(i), "time (hours)");
        %             spf2(i).XTick = 0:1/4*60*60:max(spf2(i).XLim);
        %             spf2(i).XTickLabel = spf2(i).XTick/(60*60);
        %             xlim([min(t), max(t)]);
        %         end
        %     end
        %     linkaxes_y(spf2);
        %     prepFigPresentation2(gcf)
        % 
        %     figure;
        %     sgtitle("Total Angular Momentum")
        %     for i = 3:-1:1
        %         spf3(i) = subplot(3,1,i);
        %         plot(t, Ho(:,i))
        %         grid on
        %         xlabel("time (s)")
        %         ylabel("kg*m^2/s");
        %         if i == 1
        %             text(0.9, 0.1, "ECI X", 'Units', 'normalized', 'FontSize', 16)
        %             spf3(i).XTickLabels = {};
        %         elseif i == 2
        %             text(0.9, 0.1, "ECI Y", 'Units', 'normalized', 'FontSize', 16)
        %             spf3(i).XTickLabels = {};
        %         elseif i == 3
        %             text(0.9, 0.1, "ECI Z", 'Units', 'normalized', 'FontSize', 16)
        %             xlabel(spf3(i), "time (s)");
        %         end
        %     end
        %     linkaxes_y(spf3);
        %     prepFigPresentation2(gcf)
        end

        function plotOrbitalElements(obj)
            [a, e, incl, W] = getOrbitalElements(obj);
            t = obj.simout.tout;

            figure;
            plot(t, a-a(1))
            title("Change in Semi-Major Axis")
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

        function plotDiagram(obj)
            NRac = obj.simout.NRac;
            NRbc = obj.simout.NRbc;
            NRfc = obj.simout.NRfc;
            NRgc = obj.simout.NRgc;
            NRdc = obj.simout.NRdc;
            NRec = obj.simout.NRec;

            radius = 75;
            i=1e4;

            % payload mass a
            [xa, ya, za] = sphere;
            xa = xa*radius + NRac(i, 1);
            ya = ya*radius + NRac(i, 2);
            za = za*radius + NRac(i, 3);
            
            % payload mass b
            [xb, yb, zb] = sphere;
            xb = xb*radius + NRbc(i, 1);
            yb = yb*radius + NRbc(i, 2);
            zb = zb*radius + NRbc(i, 3);

            % payload mass f
            [xf, yf, zf] = sphere;
            xf = xf*radius + NRfc(i, 1);
            yf = yf*radius + NRfc(i, 2);
            zf = zf*radius + NRfc(i, 3);

            % payload mass g
            [xg, yg, zg] = sphere;
            xg = xg*radius + NRgc(i, 1);
            yg = yg*radius + NRgc(i, 2);
            zg = zg*radius + NRgc(i, 3);

            figure;
            hold on
            axis equal
            grid on
            view(3)
            dataL = [NRac;NRbc;NRfc;NRgc;NRdc;NRec];
            axis(getAxisLimits3(dataL))
            plot3([NRdc(i,1) NRac(i,1)],[NRdc(i,2) NRac(i,2)],[NRdc(i,3) NRac(i,3)], 'r')
            plot3([NRec(i,1) NRac(i,1)],[NRec(i,2) NRac(i,2)],[NRec(i,3) NRac(i,3)], 'r')
            plot3([NRdc(i,1) NRbc(i,1)],[NRdc(i,2) NRbc(i,2)],[NRdc(i,3) NRbc(i,3)], 'r')
            plot3([NRec(i,1) NRbc(i,1)],[NRec(i,2) NRbc(i,2)],[NRec(i,3) NRbc(i,3)], 'r')
            plot3([NRdc(i,1) NRec(i,1)],[NRdc(i,2) NRec(i,2)],[NRdc(i,3) NRec(i,3)], 'k-.')
            plot3([NRdc(i,1) NRfc(i,1)],[NRdc(i,2) NRfc(i,2)],[NRdc(i,3) NRfc(i,3)], 'r')
            plot3([NRec(i,1) NRfc(i,1)],[NRec(i,2) NRfc(i,2)],[NRec(i,3) NRfc(i,3)], 'r')
            plot3([NRdc(i,1) NRgc(i,1)],[NRdc(i,2) NRgc(i,2)],[NRdc(i,3) NRgc(i,3)], 'r')
            plot3([NRec(i,1) NRgc(i,1)],[NRec(i,2) NRgc(i,2)],[NRec(i,3) NRgc(i,3)], 'r')
            plot3([NRac(i,1) NRfc(i,1)],[NRac(i,2) NRfc(i,2)],[NRac(i,3) NRfc(i,3)], 'b')
            plot3([NRac(i,1) NRgc(i,1)],[NRac(i,2) NRgc(i,2)],[NRac(i,3) NRgc(i,3)], 'b')
            plot3([NRbc(i,1) NRfc(i,1)],[NRbc(i,2) NRfc(i,2)],[NRbc(i,3) NRfc(i,3)], 'b')
            plot3([NRbc(i,1) NRgc(i,1)],[NRbc(i,2) NRgc(i,2)],[NRbc(i,3) NRgc(i,3)], 'b')
            mp_a = surf(xa, ya, za);
            mp_b = surf(xb, yb, zb);
            mp_f = surf(xf, yf, zf);
            mp_g = surf(xg, yg, zg);
            set(mp_a, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            set(mp_b, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            set(mp_f, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            set(mp_g, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
            legend(["Main Tether, Carries Current", '', '', '', "Truss", '', '', '', '', 'Safety Tether', '', '', '', 'Modules'])
            xlabel("meters")
            ylabel("meters")
            zlabel("meters")
            prepFigPresentation2(gcf)
            light('Position', [0 -1 0], 'Style', 'infinite');
        end

        function animate(obj, record)
            arguments
                obj
                record logical = false
            end

            tout = obj.simout.tout;
            tout_duration = seconds(tout);
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
            % trace_a = animatedline(axL);
            % trace_b = animatedline(axL);
            tether_ad = animatedline(axL);
            tether_ae = animatedline(axL);
            tether_bd = animatedline(axL);
            tether_be = animatedline(axL);
            % trace_f = animatedline(axL);
            % trace_g = animatedline(axL);
            tether_fd = animatedline(axL);
            tether_fe = animatedline(axL);
            tether_gd = animatedline(axL);
            tether_ge = animatedline(axL);
            truss_de = animatedline(axL);

            axR = subplot(1,2,2);
            hold on
            axis equal
            title("Spacecraft Orbit Trajectory")
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(axR,3)
            dataR = NRco;
            lims = getAxisLimits3(dataR);
            axis([lims(1:2) lims(1:2) lims(5:6)]);
            truss_center = animatedline(axR, "Color", 'r', "LineWidth", 3);
            speedUp = 20;
            [x, y, z] = sphere;
            surf(6378e3*x, 6378e3*y, 6356e3*z, 'FaceColor','b');
            hText = text(axR, 0.3, 0.85, '', 'Units','normalized');
            if record
                vidwrite = VideoWriter("spacecraft_trajectory");
                vidwrite.FrameRate = 50;
                open(vidwrite);
            end
            pause;
            for i=1:speedUp:length(NRac)
                [hour, min, ~] = hms(tout_duration(i));
                hour_disp = string(sprintf('%02d', hour));
                min_disp = string(sprintf('%02d', min));
                %sec_disp = string(sprintf('%02.0f', sec));
                hText.String = 'Time (HH:MM): +'+hour_disp+':'+min_disp;
                %addpoints(trace_a,NRac(i,1),NRac(i,2),NRac(i,3))
                %addpoints(trace_b,NRbc(i,1),NRbc(i,2),NRbc(i,3))
                addpoints(tether_ad,[NRdc(i,1) NRac(i,1)],[NRdc(i,2) NRac(i,2)],[NRdc(i,3) NRac(i,3)])
                addpoints(tether_ae,[NRec(i,1) NRac(i,1)],[NRec(i,2) NRac(i,2)],[NRec(i,3) NRac(i,3)])
                addpoints(tether_bd,[NRdc(i,1) NRbc(i,1)],[NRdc(i,2) NRbc(i,2)],[NRdc(i,3) NRbc(i,3)])
                addpoints(tether_be,[NRec(i,1) NRbc(i,1)],[NRec(i,2) NRbc(i,2)],[NRec(i,3) NRbc(i,3)])
                addpoints(truss_de,[NRdc(i,1) NRec(i,1)],[NRdc(i,2) NRec(i,2)],[NRdc(i,3) NRec(i,3)])
                %addpoints(trace_f,NRfc(i,1),NRfc(i,2),NRfc(i,3))
                %addpoints(trace_g,NRgc(i,1),NRgc(i,2),NRgc(i,3))
                addpoints(tether_fd,[NRdc(i,1) NRfc(i,1)],[NRdc(i,2) NRfc(i,2)],[NRdc(i,3) NRfc(i,3)])
                addpoints(tether_fe,[NRec(i,1) NRfc(i,1)],[NRec(i,2) NRfc(i,2)],[NRec(i,3) NRfc(i,3)])
                addpoints(tether_gd,[NRdc(i,1) NRgc(i,1)],[NRdc(i,2) NRgc(i,2)],[NRdc(i,3) NRgc(i,3)])
                addpoints(tether_ge,[NRec(i,1) NRgc(i,1)],[NRec(i,2) NRgc(i,2)],[NRec(i,3) NRgc(i,3)])
                addpoints(truss_center,NRco(i,1),NRco(i,2),NRco(i,3))
                if record && mod(i-1,speedUp*10) == 0
                    frame = getframe(gcf);
                    writeVideo(vidwrite,frame);
                end
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
            vidwrite.close();
        end
    end
end
