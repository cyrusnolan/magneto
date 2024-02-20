clc
clear
close all

% parameters
k_tether = 200;
b_tether = 20;
k_truss = 2000;
b_truss = 200;
mu = 3.986e14;
mass_payload = 3000;
mass_truss = 10;
parameters_sc = parameters(k_tether,b_tether,k_truss,b_truss,mu, ...
    mass_payload,mass_truss);

% orbit
semi_major_axis = 6878e3;
eccentricity = 0;
inclination = pi/2;
arugment_of_periapsis = 0;
lon_ascending_node = 0;
true_anomaly = pi/2;
orbit_sc = orbit(mu,semi_major_axis,eccentricity,inclination, ...
    arugment_of_periapsis,lon_ascending_node,true_anomaly);

% form spacecraft object
payload_radius = 500;
truss_length = 250;
initial_angular_velocity = 0.01;
sc = spacecraft(payload_radius,truss_length,initial_angular_velocity);

% sim spacecraft
stop_time = 300;
gravity = true;
sc = sim(sc,stop_time,gravity,parameters_sc,orbit_sc);
plot(sc)