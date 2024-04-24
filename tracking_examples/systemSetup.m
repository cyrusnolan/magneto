function [X_1dot_fun, X_ref_fun, u_ref_fun, A_lin_fun, B_lin_fun] = systemSetup(params)
%GENTRAJ generate desired tether length state and control input
%trajectories
%
    syms x y z x_1dot y_1dot z_1dot pd pe pd_1dot pe_1dot % state vars (p_1dot sym not needed)
    syms pd_2dot pe_2dot M % control inputs (p_2dot sym not needed)
    syms t 'real' % misc
    
    k = params.k;
    m = params.m;
    c = params.c;
    l0 = params.l0;
    a = params.a;
    b = params.b;
    w = params.w;
    tr = params.tr;
    g = params.g;
    phi = params.phi;
    
    % reference
    rotation = [cos(phi)  0  -sin(phi);  0  1  0;  sin(phi)  0  cos(phi)];
    Rac_ref = rotation*[a*cos(w*t); 0; b*sin(w*t)];
    Rac_1dot_ref = diff(Rac_ref, t);
    Rac_2dot_ref = diff(Rac_1dot_ref, t);
    
    % dynamics
    % position vectors
    Rac = [x; y; z];
    Rdc = -tr/2*[0; 1; 0];
    Rec = tr/2*[0; 1; 0];
    Rad = Rac - Rdc;
    Rae = Rac - Rec;
    rad = normm(Rad);
    rae = normm(Rae);
    Rad_hat = Rad/rad;
    Rae_hat = Rae/rae;
    % velocity vectors -- ignore for now
    Rac_1dot = [x_1dot; y_1dot; z_1dot];
    % mag torque direction
    M_hat = cross(Rac, Rec);
    % gravity direction
    g_hat = [0; 0; -1];
    % tension
    Td = k*(rad-l0+pd) + c*Rac_1dot.'*Rad_hat;
    Te = k*(rae-l0+pe) + c*Rac_1dot.'*Rae_hat;
    Fnet = -Td*Rad_hat - Te*Rae_hat + M*M_hat + g*g_hat;
    Rac_2dot = Fnet/m;
    % state vec
    X = [Rac; Rac_1dot; pd; pe; pd_1dot; pe_1dot];
    % input vec
    u = [pd_2dot; pe_2dot; M];
    % equations of motion
    X_1dot = [Rac_1dot; Rac_2dot; pd_1dot; pe_1dot; pd_2dot; pe_2dot];
    X_1dot_fun = matlabFunction(X_1dot, "Vars", [X; u]);
    
    % solve for reference trajectories
    Rac_2dot_eval_at_refs = subs(Rac_2dot, [Rac; Rac_1dot; pe], [Rac_ref; Rac_1dot_ref; pd]);
    acc_eqn = Rac_2dot_ref == Rac_2dot_eval_at_refs;
    acc_sol = solve(acc_eqn, pd, M);
    pd_ref = acc_sol.pd;
    M_ref = acc_sol.M;
    pd_ref = simplify(pd_ref);
    M_ref = simplify(M_ref);
    pd_1dot_ref = diff(pd_ref,t);
    pd_1dot_ref = simplify(pd_1dot_ref);
    pd_2dot_ref = diff(pd_1dot_ref,t);
    pd_2dot_ref = simplify(pd_2dot_ref);
    pe_ref = pd_ref;
    pe_1dot_ref = pd_1dot_ref;
    pe_2dot_ref = pd_2dot_ref;
    X_ref = [Rac_ref; Rac_1dot_ref; pd_ref; pe_ref; pd_1dot_ref; pe_1dot_ref];
    X_ref_fun = matlabFunction(X_ref);
    u_ref = [pd_2dot_ref; pe_2dot_ref; M_ref];
    u_ref_fun = matlabFunction(u_ref);
    
    % jacobian
    A = jacobian(X_1dot, X);
    B = jacobian(X_1dot, u);
    
    % linearize Jacobians around refs
    A_ref = subs(A, [X; M], [X_ref; M_ref]);
    A_lin_fun = matlabFunction(A_ref);
    B_ref = subs(B, X, X_ref);
    B_lin_fun = matlabFunction(B_ref);
end

