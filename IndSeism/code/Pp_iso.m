function P = Pp_iso(r, t, q, D, Ss, rho_w, g)
%%%%%
% q: injection rate, unit[m^3/s]

xi = r./sqrt(D.*t);
K = D*Ss/rho_w/g;

preFac = 1./(4*pi*K*r);
vQ_erf = q*erfc(0.5*xi);

P = preFac.*vQ_erf;
P = P/1e6;

end