clc
clear
close all
%%
% TODO
% see what body animation looks like with origin at center of mass instead of c

% BEGIN USER INPUT-----------------------------
% model = "dynamics";
tic
model = "dynamics";
mu = 3.986e14;
k_tether = 1;
b_tether = 100;
k_truss = 10;
b_truss = 100;
k_safety = 5;
b_safety = 1;
truss_length = 250;
payload_radius = 500;
mass_payload = 1500;
mass_truss = 10;
Imotor = 1;
rmotor = 1;
omega0 = 0.007;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
orbit_period = 2*pi*sqrt(semi_major_axis^3/mu);
sim_stop_time = 1*orbit_period;
% END USER INPUT--------------------------------

% initialize parameters
p = parameters(model, mu, k_tether, b_tether, k_truss, b_truss, k_safety, b_safety, ...
    truss_length, payload_radius, mass_payload, mass_truss, Imotor, rmotor, omega0);

% initialize orbit
o = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% determine initial state vector in inertial coordinates
NX = state_vec_N(p, o);

% initialize spacecraft
sc = spacecraft(p, NX);

% run sim
sc = sim(sc, sim_stop_time);
toc