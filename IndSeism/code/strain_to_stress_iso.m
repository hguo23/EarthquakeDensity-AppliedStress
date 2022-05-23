function [S_rr,S_tt,S_pp] = strain_to_stress_iso(dStrain,Pp, G, nu, biot_alpha)
%%%%%
e_rr = dStrain(1,:);
e_tt = dStrain(2,:);
e_pp = dStrain(3,:);

E = 2*G*(1+nu);
lambda = 2*G*nu/(1-2*nu);

S_rr  = E/((1+nu)*(1-2*nu))*((1-nu)*e_rr +nu*e_tt+nu*e_pp);
S_tt = (lambda+2*G)*e_tt + lambda*(e_pp+e_rr);
S_pp = (lambda+2*G)*e_pp + lambda*(e_tt+e_rr);

S_rr = S_rr/1e6-biot_alpha*Pp;
S_tt = S_tt/1e6-biot_alpha*Pp;
S_pp = S_pp/1e6-biot_alpha*Pp;

end
