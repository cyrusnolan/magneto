% Magneto Attitude Propulsion
% Initialize Simulation
% Cyrus Nolan
% May 2024

clc
clear
close all

%%
%%%%%%% BEGIN USER INPUT
% spacecraft
tether_stiffness_coeff = 10;
tether_damping_coeff = 27;
truss_stiffness_coeff = 10;
truss_damping_coeff = 10;
safety_stiffness_coeff = .1;
safety_damping_coeff = 10;
retracted_tether_length = 100;
truss_length = 250;
truss_mass = 10;
module_height = 500;
module_mass = 34.5;
spin_angular_velocity = -0.007;

% controllers
% angular momentum
% kp_amc = .0001;
kp_amc = 0.1;
zero_loc_amc = .023;
amc_on_time = 0;
% tether length
kd_tl = -0.1;
zero_loc_tl = 0.01;
tether_delta = 5;

% orbit
mu = 3.986e14;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = pi;
true_anomaly = pi/2;
orbit_period = 2*pi*sqrt(semi_major_axis^3/mu);

% misc
model = "dynamics_with_control";
% sim_stop_time = 3*orbit_period;
sim_stop_time = 2000;
%%%%%%% END USER INPUT

% initialize parameters
p = parameters(mu, tether_stiffness_coeff, tether_damping_coeff, truss_stiffness_coeff, truss_damping_coeff, safety_stiffness_coeff, safety_damping_coeff, ...
    truss_length, module_height, module_mass, truss_mass, spin_angular_velocity, tether_delta, ...
    kp_amc, zero_loc_amc, retracted_tether_length, kd_tl, zero_loc_tl, amc_on_time);

% initialize orbit
o = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% determine initial state vector in inertial coordinates
NX = state_vec_N(p, o);

% initialize spacecraft
sc = spacecraft(p, NX);

% run sim
tic
sc = sim(sc, model, sim_stop_time);
toc