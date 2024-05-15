s = tf('s');

% plant parameters
k = 10;
r_nom = 502;
l_nom = 517;
m = 1500;
alpha = r_nom / (m * l_nom);
b = 10;

% plant
Pt = 1/s^2;
Pm = -2*k*alpha / (s^2 + 2*k*alpha + s*b*alpha);
P = Pt*Pm;

control = 3
switch control
    case 1 % PID
        kp = 1;
        % Td = 1e-2;
        % Ti = 1e1;
        % ki = kp/Ti;
        % kd = kp*Td;
        kd = 1e4;
        ki = 10;
        C = (kd*s^2 + kp*s + ki) / s;
        figure;
        bode(C);

    case 2 % PI
        zero = 0.01;
        kp = 1;
        ki = zero*kp;
        C = kp * (s + ki/kp) / s;
        figure;
        bode(C);

    case 3 % PD + lead
        % PD
        zero = 0.1;
        kd = .01;
        kp = kd * zero;
        CPD = -kd * (kp/kd + s);
        % lead
        % spread = 14;
        % wL = 0.1;
        % aL = wL/sqrt(spread);
        % bL = spread*aL;
        % aL = .01;
        % bL = 1;
        % CL = bL / aL * (s + aL) / (s + bL);
        C = CPD;
        figure;
        bode(C)
        title("Controller")
        figure;
        margin(P*C)

end