clc
clear
close all

% params
p.k = 100; % spring constant
p.b = 100; % damping constant
p.ma = 1500;
p.mb = 1500;
p.mc = 10;
p.invma = 1/p.ma;
p.invmb = 1/p.mb;
p.invmc = 1/p.mc;
p.l0 = 500; % tether nominal length
p.mu = 3.986e14;

% initial conditions
ic = string_config2(p); % no gravity, no orbit
set_param("string_dynamics","StopTime","2000")
out = sim("string_dynamics.slx");

% plot
plotTraj(p,ic,out)
