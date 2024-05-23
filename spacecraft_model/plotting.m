%% Load
sim1 = sc;
%sim2 = sc2;

%% Diagram 1
plotDiagram(sc)

%% Animate 1
animate(sim1)

%% Truss 1
plotTrussLength(sim1)
 
%% Tether 1
plotTetherLength(sim1)

%% Angular Mom 1
plotAngularMom(sim1)

%% Tracking Error 1
plotTrackingError(sim1)

%% Orbital Elements 1
plotOrbitalElements(sim1)

%% Animate 2
animate(sim2)

%% Truss 2
plotTrussLength(sim2)

%% Tether 2
plotTetherLength(sim2)

%% Angular Mom 2
plotAngularMom(sim2)

%% Orbital Elements 2
plotOrbitalElements(sim2)
