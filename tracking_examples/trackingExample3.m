% simplified magneto 3
% track a reference x^2/a^2 + y^2/b^2 = 1
% same as simplified magneto 2, but gravity is now included

clc
clear
close all
%%
% params
a = 3; p.a = a;
b = 2; p.b = b;
m = 1; p.m = m;
g = 9.8; p.g = g;
w = 0.5;

% symbolic system
[Xdot_fun, A_lin_fun, B_lin_fun, Xd_fun, ud_fun] = sytemSetup(p);

% simulate
tspan = 0:0.1:20;
X0 = [a*cos(0);  b*sin(0);  -a*w*sin(0);  b*w*cos(0)];
[tout, yout] = ode45(@(t, x) odefun(t, x, Xdot_fun, A_lin_fun(t,w), B_lin_fun(t,w), ...
    Xd_fun(t,w), ud_fun(t,w)), tspan, X0);

% plot
x = yout(:,1);
y = yout(:,2);
plot(x,y)
axis equal


function [Xdot_fun, A_lin_fun, B_lin_fun, Xd_fun, ud_fun] = sytemSetup(p)
    % vars
    syms x xdot y ydot F T l t w 'real'

    % params
    a = p.a;
    b = p.b;
    m = p.m;
    g = p.g;
    
    % O coordinate system
    o1 = [x/l; y/l];
    o2 = [-y/l; x/l];
    
    % state vec
    r = [x; y];
    rdot = [xdot; ydot];
    X = [r; rdot];
    
    % input vec
    u = [F; T];
    
    % equations of motion
    rddot = -F/m*o1 + T/m*o2 - g*[0; 1];
    Xdot = [rdot; rddot];
    Xdot_fun = matlabFunction(Xdot, "Vars", {F, T ,l, x, xdot, y, ydot});
    
    % jacobian
    A = jacobian(Xdot,X);
    B = jacobian(Xdot,u);

    % reference state trajectory
    xd = a*cos(w*t);
    yd = b*sin(w*t);
    ld = sqrt(xd^2 + yd^2);
    xdotd = -a*w*sin(w*t);
    ydotd = b*w*cos(w*t);
    Xd = [xd; yd; xdotd; ydotd];
    Xd_fun = matlabFunction(Xd, "Vars", {t,w});

    % reference control trajectory
    xddotd = -a*w^2*cos(w*t);
    yddotd = -b*w^2*sin(w*t);
    rddotd = subs(rddot, [x, y, l], [xd, yd, ld]);
    syseqn = [xddotd; yddotd] == rddotd;
    [AA,BB] = equationsToMatrix(syseqn,[F,T]);
    ud = linsolve(AA,BB);
    ud_fun = matlabFunction(ud, "Vars", {t,w});

    % evaluate jacobians at xd ud
    Fd = ud(1);
    Td = ud(2);

    A_lin = subs(A, [F, T, l], [Fd, Td, ld]);
    A_lin_fun = matlabFunction(A_lin, "Vars", {t,w});
    B_lin = subs(B, [l, x, y], [ld, xd, yd]);
    B_lin_fun = matlabFunction(B_lin, "Vars", {t,w});
end

function Xdot = odefun(~, X, Xdot_fun, A, B, Xd, ud)
    % X = [x  y  xdot  ydot]'
    % u = [F  T];
    
    % state vec
    x = X(1);
    y = X(2);
    xdot = X(3);
    ydot = X(4);

    % solve lqr
    Q = eye(4);
    R = eye(2);
    [K,~,~] = lqr(A,B,Q,R);

    % control law
    u = ud - K*(X-Xd);

    % progate state
    l = sqrt(x^2 + y^2);
    F = u(1);
    T = u(2);
    Xdot = Xdot_fun(F, T, l, x, xdot, y, ydot);
end