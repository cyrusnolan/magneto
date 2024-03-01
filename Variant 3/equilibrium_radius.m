function re = equilibrium_radius(parameters, w)
% equilibrium calculation
k = paramters.k;
l0 = parameters.l0;

syms r
Fa = 2*
FC = -p.ma*(w0*r)^2/r;
eqn = Fa - FC == 0;
re = double(solve(eqn));
end

