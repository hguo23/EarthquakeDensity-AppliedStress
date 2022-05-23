function P = Pp_axi(r, t, q, D, Ss, b, rho_w, g)
%%%%%
% q: injection rate, unit[m^3/s]

T = D*Ss*b;

u = r.^2./(4*D*t');
P = q*rho_w*g/(4*pi*T)*expint(u);

P = P/1e6;
end
