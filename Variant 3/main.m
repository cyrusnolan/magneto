clc
clear
close all

% initialize at equilibrium position and add expectation equilibrium to
% plot() figures 1 and 2

% parameters
mu = 3.986e14;
k_tether = 10;
b_tether = 100;
k_truss = 300;
b_truss = 200;
truss_length = 250;
tether_length = sqrt(500^2+truss_length^2);
mass_payload = 3000;
mass_truss = 10;
equilibrium_angular_velocity = 0.01;
v3_parameters = parameters(mu, k_tether, b_tether, k_truss, b_truss, truss_length, tether_length, mass_payload, ...
    mass_truss, equilibrium_angular_velocity);

% orbit
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
v3_orbit = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% form spacecraft
sc = spacecraft(v3_parameters);
sc = equilibrium(sc);


% % sim spacecraft
% stop_time = 7000;
% gravity = true;
% sc = sim(sc, stop_time, gravity, v3_orbit);
% % plot(sc)
% animate(sc)