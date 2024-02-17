function re = equilibrium_radius(w0,p)
% equilibrium calculation
syms r
FT = -p.k*(r-p.l0);
FC = -p.ma*(w0*r)^2/r;
eqn = FT - FC == 0;
re = double(solve(eqn));
end

