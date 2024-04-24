% simplified magneto 1
% track a reference ellipse

clc
clear
close all
%%
params.a = 10;
params.b = 8;
params.m = 10;
params.w = 0.5;
params.k = 10;
params.c = 1;
params.l0 = 11;
params.tr = 4;
params.g = 9.8;
params.phi = deg2rad(45);

[X_1dot_fun, X_ref_fun, u_ref_fun, A_lin_fun, B_lin_fun] = systemSetup(params);

tspan = 0:0.01:20;
X0 = X_ref_fun(0);

[tout, yout] = ode45(@(t, X) ...
    odefun(t, X, X_1dot_fun, X_ref_fun(t), u_ref_fun(t), A_lin_fun(t), B_lin_fun(t)), ...
    tspan, X0);

% eval things
x = yout(:,1);
y = yout(:,2);
z = yout(:,3);
r = (x.^2 + y.^2 + z.^2).^(1/2);

% plot things
figure;
plot3(x, y, z);
axis equal;
title("trajectory")
xlabel("x (m)")
ylabel("y (m)")
zlabel("z (m)")

figure;
plot(tout, r)
title("tether length")
xlabel("time (s)")
ylabel("length (m)")

function Xdot = odefun(~, X, X_1dot_fun, X_ref, u_ref, A, B)
    % solve lqr
    Q = eye(10);
    R = eye(3);
    [K,~,~] = lqr(A,B,Q,R);
    
    % control law
    u = u_ref - K*(X-X_ref);
    
    % propagate state
    args = num2cell([X; u]);
    Xdot = X_1dot_fun(args{:});
end