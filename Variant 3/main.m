clc
clear
close all

% TODO
% change some/all angular mom plots to normalized by initial value
% remove parameter config - deprecated

% BEGIN USER INPUT-----------------------------
config = "octahedron"; % config = {"v1", "octahedron"}
sim_stop_time = 86400;
mu = 3.986e14;
k_tether = 10;
b_tether = 100;
k_truss = 10;
b_truss = 100;
k_safety = 5;
b_safety = 1;
truss_length = 250;
payload_radius = 500;
mass_payload = 1500;
mass_truss = 10;
omega0 = -0.01;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
% END USER INPUT--------------------------------

% initialize parameters
p = parameters(config, mu, k_tether, b_tether, k_truss, b_truss, k_safety, b_safety, ...
    truss_length, payload_radius, mass_payload, mass_truss, omega0);

% initialize orbit
orbit = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% determine initial state vector in inertial coordinates
NX = state_vec_N(p, orbit);

% initialize spacecraft
sc = spacecraft(p, NX);

% run sim
sc = sim(sc, sim_stop_time);