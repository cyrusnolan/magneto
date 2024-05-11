clc
clear
close all

%%
% BEGIN USER INPUT-----------------------------
% model = "dynamics";
% next step: plot phi on tether length plot to see if tether length long
% where I want it.
tic
model = "dynamics_with_control";
mu = 3.986e14;
p0 = 100;
kp_M = .1;
zero_loc_M = .023;
k_tether = 10;
b_tether = 10;
k_truss = 10;
b_truss = 10;
k_safety = .1;
b_safety = 10;
truss_length = 250;
payload_radius = 500;
tether_delta = 4;
mass_payload = 1500;
mass_truss = 10;
omega0 = -0.007;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = pi;
true_anomaly = pi/2;
orbit_period = 2*pi*sqrt(semi_major_axis^3/mu);
% sim_stop_time = 1*orbit_period;
sim_stop_time = 3000;
% END USER INPUT--------------------------------

% initialize parameters
p = parameters(model, mu, k_tether, b_tether, k_truss, b_truss, k_safety, b_safety, ...
    truss_length, payload_radius, mass_payload, mass_truss, omega0, tether_delta, ...
    kp_M, zero_loc_M, p0);

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