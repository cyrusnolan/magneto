% simplified magneto 4
% track rotated ellipse
% same as simplified magneto 3, but now the ellipse is tilted, the dynamics include damping
% *note names of input forces have changed from example 3

clc
clear
close all
%%
% params
p.a = 550;
p.b = 450;
p.ma = 1;
p.g = 9.8;
p.k = 100;
p.c = 10;
p.w = -0.5;
p.phi = deg2rad(45);

% symbolic system
[Xdot_fun, Ad_fun, Bd_fun, Xd_fun, ud_fun] = sytemSetup(p);

% simulate
tspan = 0:0.1:20;
X0 = Xd_fun(0);
[tout, yout] = ode45(@(t, x) odefun(t, x, Xdot_fun, Ad_fun(t), Bd_fun(t), Xd_fun(t), ud_fun(t)), ...
                                tspan, X0);

% plot desired inputs
for i = length(tspan):-1:1
    ud_eval(1:2,i) = ud_fun(tspan(i));
end
figure;
subplot(1,2,1)
plot(tspan, ud_eval(1,:))
title("Tud")
xlabel("time (s)")
ylabel("Tud (N)")

subplot(1,2,2)
plot(tspan, ud_eval(2,:))
title("Md")
xlabel("time (s)")
ylabel("Md (N)")

% plot output
x = yout(:,1);
y = yout(:,2);
figure;
plot(x,y)
title("Output trajectory")
xlabel("x (m)")
ylabel("y (m)")
axis equal

function [Xdot_fun, Ad_fun, Bd_fun, Xd_fun, ud_fun] = sytemSetup(p)
    % vars
    syms x xdot y ydot l0 l0dot 'real' % state variables
    syms Tu M 'real' % control variables
    syms t l 'real' % misc
    
    % params
    a = p.a;
    b = p.b;
    ma = p.ma;
    g = p.g;
    k = p.k;
    c = p.c;
    w = p.w;
    phi = p.phi;
        
    % coordinate systems
    n2 = [0; 1];
    o1 = [x/l; y/l];
    o2 = [-y/l; x/l];
        
    % state vec
    r = [x; y];
    rdot = [xdot; ydot];
    X = [r; rdot; l0; l0dot];
    
    % input vec
    u = [Tu; M];
    
    % equations of motion
    T = k*(l - l0) + c*dot(rdot,o1);
    rddot = 1/ma*(-T*o1 + M*o2 - g*n2);
    l0ddot = T-Tu;
    Xdot = [rdot; rddot; l0dot; l0ddot];
    Xdot_fun = matlabFunction(Xdot, "Vars", {x, xdot, y, ydot, l0, l0dot, l, Tu, M});
    
    % jacobian
    A = jacobian(Xdot,X);
    B = jacobian(Xdot,u);
    
    % reference trajectories
    % desired mass position, velocity, acceleration
    NQB = [cos(phi)  -sin(phi);  sin(phi)  cos(phi)];
    Brd = [a*cos(w*t);  b*sin(w*t)];
    rd = NQB*Brd;
    xd = rd(1);
    yd = rd(2);
    rdotd = diff(rd, t);
    xdotd = rdotd(1);
    ydotd = rdotd(2);
    rddotd = diff(rdotd, t);
    
    % mass acceleration must satisfy equations of motion
    % generates desired tether length (l0d) and desired magnetic force (Md)
    ld = sqrt(xd^2 + yd^2);
    r_dynamics_d = subs(rddot, [x, y, xdot, ydot, l], [xd, yd, xdotd, ydotd, ld]);
    syseqn = rddotd == r_dynamics_d;
    [AA, BB] = equationsToMatrix(syseqn, [l0, M]);
    sol = linsolve(AA, BB);
    l0d = sol(1);
    Md = sol(2);
    
    % tether length acceleration must satisfy equation of motion
    % generates desired force on tether by motor (Tud)
    l0dotd = diff(l0d,t);
    l0ddotd = diff(l0dotd,t);
    l0_dynamics_d = subs(l0ddot, [x, y, xdot, ydot, l, l0], [xd, yd, xdotd, ydotd, ld, l0d]);
    eqn = l0ddotd == l0_dynamics_d;
    Tud = solve(eqn,Tu);
    
    % package up desired trajectories
    Xd = [xd; yd; xdotd; ydotd; l0d; l0dotd];
    Xd_fun = matlabFunction(Xd, "Vars", t);
    ud = [Tud; Md];
    ud_fun = matlabFunction(ud, "Vars", t);
    
    % evaluate jacobians at xd, ud
    Ad = subs(A, [x, y, xdot, ydot, l, l0, M], [xd, yd, xdotd, ydotd, ld, l0d, Md]);
    Ad_fun = matlabFunction(Ad, "Vars", t);
    Bd = subs(B, [x, y, l], [xd, yd, ld]);
    Bd_fun = matlabFunction(Bd, "Vars", t);
end

function Xdot = odefun(~, X, Xdot_fun, A, B, Xd, ud)
    % X = [x  y  xdot  ydot  l0  l0dot]'
    % u = [Tu  M];
    
    % state vec
    x = X(1);
    y = X(2);
    xdot = X(3);
    ydot = X(4);
    l0 = X(5);
    l0dot = X(6);
    
    % solve lqr
    Q = eye(6);
    R = eye(2);
    [K,~,~] = lqr(A,B,Q,R);
    
    % control law
    u = ud - K*(X-Xd);
    
    % propagate state
    l = sqrt(x^2 + y^2);
    Tu = u(1);
    M = u(2);
    Xdot = Xdot_fun(x, xdot, y, ydot, l0, l0dot, l, Tu, M);
      % arg order: {x, xdot, y, ydot, l0, l0dot, l, Tu, M}
end