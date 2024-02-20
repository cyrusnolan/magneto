clc
clear
close all
% finish spacecraft class
% simulate dynamics model asap and start debugging

p1 = parameters(1000,20,3.986e14,3000,20);
o1 = orbit(3.986e14,6878e3,0,pi/2,0,0,pi/2);
sc = spacecraft(500,250,0.05);
sc = sim(sc,300,true,p1,o1);
plot(sc)