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