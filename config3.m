function ic = config3()
% 3d orbit with gravity
% a = 6878km, e = 0, i = 90 deg, Raan = 180 deg
% --- set these ---
ic.gravity = true;
R1 = deg2rad(0); % zxz rotation angles
R2 = deg2rad(90);
R3 = 0;
BWc0 = [0;0;.01]; % w_/c = w_a/c = w_b/c
BRac0 = 500*[1;0;0]; % r_a/c
BRbc0 = 500*[-1;0;0]; % r_b/c
NRco0 = 6878e3*[0;0;1]; % r_c/o
ic.NRco0 = NRco0;
NVco0 = sqrt(3.986e14/norm(NRco0))*[1;0;0]; % v_c/o
ic.NVco0 = NVco0;
% -----------------
NQF = [cos(R1) -sin(R1) 0;sin(R1) cos(R1) 0;0 0 1];
FQG = [1 0 0;0 cos(R2) -sin(R2);0 sin(R2) cos(R2)];
GQB = [cos(R3) -sin(R3) 0;sin(R3) cos(R3) 0;0 0 1];
NQB = NQF*FQG*GQB;
BVac0 = cross(BWc0,BRac0);
BVbc0 = cross(BWc0,BRbc0);
NRac0 = NQB*BRac0;
NVac0 = NQB*BVac0;
NRbc0 = NQB*BRbc0;
NVbc0 = NQB*BVbc0;
ic.NRao0 = NRac0 + NRco0; % ma initial position
ic.NVao0 = NVac0 + NVco0; % ma initial velocity
ic.NRbo0 = NRbc0 + NRco0; % mb initial position
ic.NVbo0 = NVbc0 + NVco0; % mb initial velocity
end

