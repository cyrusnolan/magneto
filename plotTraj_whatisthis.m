function [] = plotTraj(Rac,Rbc,Rco)
xa = Rac(:,1);
ya = Rac(:,2);
za = Rac(:,3);

xb = Rbc(:,1);
yb = Rbc(:,2);
zb = Rbc(:,3);

xc = Rco(:,1);
yc = Rco(:,2);
zc = Rco(:,3);

xrange = max([xa;xb;xc]) - min([xa;xb;xc]);
yrange = max([ya;yb;yc]) - min([ya;yb;yc]);
zrange = max([za;zb;zc]) - min([za;zb;zc]);
xmin = min([xa;xb]) - 0.1*xrange;
xmax = max([xa;xb]) + 0.1*xrange;
ymin = min([ya;yb]) - 0.1*yrange;
ymax = max([ya;yb]) + 0.1*yrange;
zmin = min([za;zb]) - 0.1*zrange;
zmax = max([za;zb]) + 0.1*zrange;

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
xlim([xmin xmax])
ylim([ymin ymax])
zlim([zmin zmax])
xlabel("ECI X")
ylabel("ECI Y")
zlabel("ECI Z")
plot3(axL,xa,ya,za)
plot3(axL,xb,yb,zb)
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
plot3(axR,xc,yc,zc)
end

