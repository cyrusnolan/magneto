% simplified magneto 1
% track a reference ellipse

clc
clear
close all
%%
% params
params.a = 10;
params.b = 5;
params.m = 1;
params.w = 0.5;
params.k = 10;
params.c = 1;
params.l0 = 10;

% system as a function of time
[X_1dot_fun, A_lin_fun, B_lin_fun, X_ref_fun, u_ref_fun] = sytemSetup(params);

% simulate
tspan = 0:0.1:20;
X0 = X_ref_fun(0);
[tout, yout] = ode45(@(t, X) odefun(t, X, X_1dot_fun, A_lin_fun(t), B_lin_fun(), ...
    X_ref_fun(t), u_ref_fun(t)), tspan, X0);

% plot
x = yout(:,1);
y = yout(:,2);
z = yout(:,3);
plot3(x, y, z);


function [X_1dot_fun, A_lin_fun, B_lin_fun, X_ref_fun, u_ref_fun] = sytemSetup(params)
    % vars
    syms x y z x_1dot y_1dot z_1dot p p_1dot % state vars
    syms p_2dot % control input
    syms t 'real'
    
    % params
    a = params.a;
    b = params.b;
    m = params.m;
    w = params.w;
    k = params.k;
    c = params.c;
    l0 = params.l0;
    
    % state vec
    Rac = [x; y; z];
    Rac_1dot = [x_1dot; y_1dot; z_1dot];
    X = [Rac; Rac_1dot; p; p_1dot];
    
    % input vec
    u = p_2dot;
    
    % equations of motion
    Rac_hat = Rac/normm(Rac);
    Rac_2dot = 1/m*(-(k*(l - l0 - p) + c*dott(Rac_1dot, Rac_hat))*Rac_hat);
    X_1dot = [Rac_1dot; Rac_2dot; p_1dot; p_2dot];
    X_1dot = subs(X_1dot, l, normm(Rac));
    X_1dot_fun = matlabFunction(X_1dot, "Vars", [X; u]);
    
    % jacobian
    A = jacobian(X_1dot, X);
    B = jacobian(X_1dot, u);
    
    % reference trajectory
    P_ref = tether_length(params);
    p_ref = P_ref(1);
    p_1dot_ref = P_ref(2);
    p_2dot_ref = P_ref(3);
    Rac_ref = [a*cos(w*t); 0; b*sin(w*t)];
    Rac_1dot_ref = diff(Rac_ref, t);
    l_ref = normm(Rac_ref);
    
    % evaluate jacobians at xd ud
    A_lin = subs(A, [Rac; Rac_1dot; p; l], [Rac_ref; Rac_1dot_ref; p_ref; l_ref]);
    A_lin_fun = matlabFunction(A_lin, "Vars", t);
    B_lin = B;
    B_lin_fun = matlabFunction(B_lin);
    
    % package trajectories for return
    % state trajectories
    X_ref = [Rac_ref; Rac_1dot_ref; p_ref; p_1dot_ref];
    X_ref_fun = matlabFunction(X_ref);
    % control input trajectory
    u_ref = p_2dot_ref;
    u_ref_fun = matlabFunction(u_ref);
end

function Xdot = odefun(~, X, X_1dot_fun, A, B, X_ref, u_ref)
    % solve lqr
    Q = eye(8);
    R = eye(1);
    [K,~,~] = lqr(A,B,Q,R);

    % control law
    u = u_ref - K*(X-X_ref);

    % propagate state
    Xdot = X_1dot_fun(X, u);
end