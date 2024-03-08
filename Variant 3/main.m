clc
clear
close all

% TODO
% add tension forces to mass d and e in dynamics_v2 simulink model
% decrease function tolerance on equilibrium calculation fsolve

% BEGIN USER INPUT
mu = 3.986e14;
k_tether = 10;
b_tether = 100;
k_truss = 300;
b_truss = 300;
truss_length = 250;
payload_radius = 500;
num_payload = 4; % 2 payload mass or 4 payload mass configuration
mass_payload = num_payload/2*3000;
mass_truss = 10;
initial_angular_velocity = 0.01;
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
% END USER INPUT

% initialize parameters
v1_parameters = parameters(mu, k_tether, b_tether, k_truss, b_truss, truss_length, payload_radius, mass_payload, ...
    num_payload, mass_truss, initial_angular_velocity);

% determine initial state vector
v1_state = state_vec(v1_parameters);

% initialize spacecraft
v1 = spacecraft(v1_parameters, v1_state);

% initialize orbit
v1_orbit = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% sim spacecraft
stop_time = 8000;
v1 = sim(v1, stop_time, v1_orbit);
plot(v1)
% animate(v1)