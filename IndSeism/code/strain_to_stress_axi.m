function [S_rr,S_tt,S_zz] = strain_to_stress_axi(dStrain,Pp, G, nu, biot_alpha)
%%%%%
e_rr = dStrain(1,:);
e_tt = dStrain(2,:);
e_kk = dStrain(3,:);

S_rr = 2*G*e_rr+2*G*nu./(1-2*nu)*e_kk;
S_tt = 2*G*e_tt+2*G*nu./(1-2*nu)*e_kk;

S_rr = S_rr/1e6-biot_alpha*Pp;
S_tt = S_tt/1e6-biot_alpha*Pp;
S_zz = nu*(S_rr+S_tt)-(1-2*nu)*biot_alpha*Pp;
end


        
