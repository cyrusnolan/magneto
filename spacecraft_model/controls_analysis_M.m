s = tf('s');

% plant parameters
k = 10;
r_nom = 502;
l_nom = 517;
wn = 0.75;
m = 2*k/(wn^2)*r_nom/l_nom;
zeta = .5;  
b = 2*zeta*wn*l_nom*m/r_nom;
alpha = r_nom / (m * l_nom);

% plant
Pt = 1/s^2;
Pm = -2*k*alpha / (s^2 + 2*k*alpha + s*b*alpha);
P = Pt*Pm;

control = 3;
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
        zero = .01;
        kd = 0.1;
        kp = kd * zero;
        CPD = -kd * (kp/kd + s);
        C = CPD;
        % C = -1;
        figure;
        margin(P*C)
        prepBodePresentation(gcf)
end