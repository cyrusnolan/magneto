function [] = animateTraj2(p,ic,out)
% TODO

% read in coords
Rao = out.Rao;
Rbo = out.Rbo;
Rco = out.Rco;
Rac = zeros(size(Rao));
Rbc = zeros(size(Rao));
rac = zeros(length(Rao),1);
rbc = zeros(length(Rbo),1);
for i = 1:length(Rao)
    Rac(i,:) = Rao(i,:) - Rco(i,:);
    Rbc(i,:) = Rbo(i,:) - Rco(i,:);
    rac(i) = norm(Rac(i,:));
    rbc(i) = norm(Rbc(i,:));
end
xao = Rao(:,1);
yao = Rao(:,2);
zao = Rao(:,3);
xbo = Rbo(:,1);
ybo = Rbo(:,2);
zbo = Rbo(:,3);
xco = Rco(:,1);
yco = Rco(:,2);
zco = Rco(:,3);
xac = Rac(:,1);
yac = Rac(:,2);
zac = Rac(:,3);
xbc = Rbc(:,1);
ybc = Rbc(:,2);
zbc = Rbc(:,3);

% find maxes and mins of Rco (xco, yco, zco)
xorange = max(xco) - min(xco);
yorange = max(yco) - min(yco);
zorange = max(zco) - min(zco);
xomin = min(xco) - 0.1*xorange;
xomax = max(xco) + 0.1*xorange;
yomin = min(yco) - 0.1*yorange;
yomax = max(yco) + 0.1*yorange;
zomin = min(zco) - 0.1*zorange;
zomax = max(zco) + 0.1*zorange;

% find maxes and mins of Rac and Rbc ([xac;xbc], [yac;ybc], [zac;zbc])
xcrange = max([xac;xbc]) - min([xac;xbc]);
ycrange = max([yac;ybc]) - min([yac;ybc]);
zcrange = max([zac;zbc]) - min([zac;zbc]);
xcmin = min([xac;xbc]) - 0.1*xcrange;
xcmax = max([xac;xbc]) + 0.1*xcrange;
ycmin = min([yac;ybc]) - 0.1*ycrange;
ycmax = max([yac;ybc]) + 0.1*ycrange;
zcmin = min([zac;zbc]) - 0.1*zcrange;
zcmax = max([zac;zbc]) + 0.1*zcrange;

% plot 1: 3d trajectory
axL = subplot(2,2,[1,3]);
grid on
axis equal
hold on
[oranges,idx] = sort([xorange yorange zorange]);
if idx(1) == 1
    xlim([min(xco)-0.1*oranges(2), max(xco)+0.1*oranges(2)])
    ylim([yomin yomax])
    zlim([zomin zomax])
elseif idx(1) == 2
    xlim([xomin xomax])
    ylim([min(yco)-0.1*oranges(2), max(yco)+0.1*oranges(2)])
    zlim([zomin zomax])
elseif idx(1) == 3
    xlim([xomin xomax])
    ylim([yomin yomax])
    zlim([min(zco)-0.1*oranges(2), max(zco)+0.1*oranges(2)])
end
view(3)
xlabel("ECI X (m)")
ylabel("ECI Y (m)")
zlabel("ECI Z (m)")
title("Trajectories relative to O")
% legendL = ["Mass a","Mass b","Mass c"];
% legend(legendL,"Location","northwest")

% plot 2: 2d trajectory
axRT = subplot(2,2,2);
grid on
axis equal
hold on
view(ic.NWc0)
xlabel("ECI X (m)")
ylabel("ECI Y (m)")
zlabel("ECI Z (m)")
title("Trajectory of mass a and b relative to c")
% legendRT = ["Mass a","Mass b"];
% legend(legendRT)

% plot 3: magnitude of payload radius
axRB = subplot(2,2,4);
grid on
axis equal
re = equilibrium_radius(norm(ic.BWc0),p);
ylim(axRB,[re-p.l0-1,re-p.l0+1])
ylim(axRB,[re-p.l0-1,re-p.l0+1])
title("Magnitude of position relative to spacecraft center")
% legendRB = ["Equilibrium","Mass a","Mass b"];
% legend(legendRB)

for i = 1:length(xao)
    tethera = plot3(axRT,[0 xac(i)],[0 yac(i)],[0 zac(i)],'k',LineWidth=1.2);
    tetherb = plot3(axRT,[0 xbc(i)],[0 ybc(i)],[0 zbc(i)],'r',LineWidth=1.2);
    plot3(axL,xco(1:i),yco(1:i),zco(1:i),'k');
    plot(axRB,out.tout(1:i),rac(1:i)-p.l0,"r")
    plot(axRB,out.tout(1:i),rbc(1:i)-p.l0,"b-.")
    pause(.001)
    delete(tethera);
    delete(tetherb);
end
