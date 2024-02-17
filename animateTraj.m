function [] = animateTraj(Rac,Rbc,Rco)
xa = Rac(:,1);
ya = Rac(:,2);
za = Rac(:,3);

xb = Rbc(:,1);
yb = Rbc(:,2);
zb = Rbc(:,3);

xc = Rco(:,1);
yc = Rco(:,2);
zc = Rco(:,3);

xLrange = max([xa;xb]) - min([xa;xb]);
yLrange = max([ya;yb]) - min([ya;yb]);
zLrange = max([za;zb]) - min([za;zb]);
xLmin = min([xa;xb]) - 0.1*xLrange;
xLmax = max([xa;xb]) + 0.1*xLrange;
yLmin = min([ya;yb]) - 0.1*yLrange;
yLmax = max([ya;yb]) + 0.1*yLrange;
zLmin = min([za;zb]) - 0.1*zLrange;
zLmax = max([za;zb]) + 0.1*zLrange;
xRrange = max(xc) - min(xc);
yRrange = max(yc) - min(yc);
zRrange = max(zc) - min(zc);
xRmin = min(xc) - 0.1*xRrange;
xRmax = max(xc) + 0.1*xRrange;
yRmin = min(yc) - 0.1*yRrange;
yRmax = max(yc) + 0.1*yRrange;
zRmin = min(zc) - 0.1*zRrange;
zRmax = max(zc) + 0.1*zRrange;

axL = subplot(1,2,1);
view(3)
grid on
axis equal
hold on
xlim([xLmin xLmax])
ylim([yLmin yLmax])
zlim([zLmin zLmax])
xlabel("ECI X")
ylabel("ECI Y")
zlabel("ECI Z")
axR = subplot(1,2,2);
view(3)
grid on
axis equal
hold on
xlim([xRmin xRmax])
ylim([yRmin yRmax])
zlim([zRmin zRmax])
xlabel("ECI X")
ylabel("ECI Y")
zlabel("ECI Z")

for i = 1:length(xa)
    tethera = plot3(axL,[0 xa(i)],[0 ya(i)],[0 za(i)],'k',LineWidth=1.2);
    tetherb = plot3(axL,[0 xb(i)],[0 yb(i)],[0 zb(i)],'r',LineWidth=1.2);
    if i > 1
        tracea = plot3(axL,xa(i-1:i),ya(i-1:i),za(i-1:i),'k');
        traceb = plot3(axL,xb(i-1:i),yb(i-1:i),zb(i-1:i),'r');
        tracec = plot3(axR,xc(i-1:i),yc(i-1:i),zc(i-1:i),'k');
    end
    pause(.001)
    delete(tethera);
    delete(tetherb);
end
end
