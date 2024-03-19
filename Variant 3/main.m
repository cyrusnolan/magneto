clc
clear
close all

% TODO
% 

% BEGIN USER INPUT-----------------------------
sim_stop_time = 80000;
model = 'octahedron'; % {"v1", "octahedron"}
mu = 3.986e14;
k_tether = 10;
b_tether = 100;
k_truss = 10;
b_truss = 100;
truss_length = 250;
payload_radius = 500;
mass_payload = 1500;
mass_truss = 10;
omega0 = 0.01;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
% END USER INPUT--------------------------------

% initialize parameters
p = parameters(model, mu, k_tether, b_tether, k_truss, b_truss, truss_length, payload_radius, ...
    mass_payload, mass_truss, omega0);

% determine initial state vector in body coordinates
BX = state_vec_B(p);

% initialize orbit
orbit = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% determine initial state vector in inertial coordinates
NX = state_vec_N(BX, orbit);

% initialize spacecraft
sc = spacecraft(p, NX);

% run sim
sc = sim(sc, sim_stop_time);