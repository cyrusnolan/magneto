% simplified magneto 1
% track a reference x^2/a^2 + y^2/b^2 = 1
% one control input - radial force applied directly to the mass

clc
clear
close all
%%
% params
syms t x xdot y ydot u l 'real'
a = 3; p.a = a;
b = 1; p.b = b;
w = .5; p.w = w;
m = 1; p.m = m;

% symbolic A and B
r = [x; y];
rdot = [xdot; ydot];
rddot = -u/m*[x/l; y/l];
X_s = [r; rdot];
Xdot_s = [rdot; rddot];
A_s = jacobian(Xdot_s,X_s);
B_s = jacobian(Xdot_s,u);

% evaulate at xd, ud
xd = [a*cos(w*t)  b*sin(w*t)];
ud = m*w^2*sqrt(a^2*cos(w*t)^2 + b^2*sin(w*t)^2);
A_lin_s = subs(A_s, [x  y  u  l], [xd  ud  norm(xd)]);
A = simplify(A_lin_s);
B_lin_s = subs(B_s, [x  y  l], [xd  norm(xd)]);
B_lin = matlabFunction(B_lin_s, "Vars", {t});

% analyze stabalizability with PBH test
% linearized A is LTI
ti = 0;
tf = 4*pi/w;
dt = pi/32;
tspan = ti:dt:tf;
L = eig(A);

for i = length(tspan):-1:1
    B = B_lin(tspan(i));
    PBH1(i) = rank([A-eye(4)*L(1)  B]);
end

plot(tspan, PBH1) % need another control input