function ic = string_config2(p)
% yes gravity, yes orbit
% SI units
% N indicates ECI coordinates
% B indicates body-fixed coordiates

% --- set these ---
ic.gravity = true;
R1 = deg2rad(90); % zxz rotation angles
R2 = deg2rad(90);
R3 = deg2rad(0);
wc0 = .05; % magnitude of angular velocity about c
BWc0hat = [0;0;1]; % direction of angular velocity about c
rac0 = equilibrium_radius(wc0,p); % magnitude of radius a/c
rbc0 = equilibrium_radius(wc0,p); % magnitude of radius b/c
rco0 = 6878e3; % magnitude of radius c/o
NRco0hat = [1;0;0]; % direction of radius c/o
vco0 = sqrt(p.mu/rco0); % magnitude of velocity c/o
NVco0hat = [0;1;0]; % direction of velocity c/o
% -----------------

ic.BWc0 = wc0*BWc0hat; % w_/c = w_a/c = w_b/c
BRac0 = rac0*[1;0;0]; % r_a/c
BRbc0 = rbc0*[-1;0;0]; % r_b/c
ic.NRco0 = rco0*NRco0hat; % r_c/o
ic.NVco0 = vco0*NVco0hat; % v_c/o
NQF = [cos(R1) -sin(R1) 0;sin(R1) cos(R1) 0;0 0 1];
FQG = [1 0 0;0 cos(R2) -sin(R2);0 sin(R2) cos(R2)];
GQB = [cos(R3) -sin(R3) 0;sin(R3) cos(R3) 0;0 0 1];
NQB = NQF*FQG*GQB;
BVac0 = cross(ic.BWc0,BRac0);
BVbc0 = cross(ic.BWc0,BRbc0);
NRac0 = NQB*BRac0;
NVac0 = NQB*BVac0;
NRbc0 = NQB*BRbc0;
NVbc0 = NQB*BVbc0;
ic.NWc0 = NQB*ic.BWc0;
ic.NRao0 = NRac0 + ic.NRco0; % ma initial position
ic.NVao0 = NVac0 + ic.NVco0; % ma initial velocity
ic.NRbo0 = NRbc0 + ic.NRco0; % mb initial position
ic.NVbo0 = NVbc0 + ic.NVco0; % mb initial velocity

end

