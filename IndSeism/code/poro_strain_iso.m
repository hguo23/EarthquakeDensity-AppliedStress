function dStrain = poro_strain_iso(r, t, q, D, alpha, G, lambda, lambda_u,varargin)
%%%%%
r = reshape(r,1,length(r));
preFac = (lambda_u-lambda)/(8*pi*D*alpha*(lambda_u+2*G));

%%
dr = 1e-5;
rvec = [r-dr;r;r+dr];
disp = nan(3,length(r));
for i = 1:length(r)
    xi = rvec(:,i)/sqrt(D*t);   
    g1 = erf(0.5*xi) - (xi/sqrt(pi)).*exp(-.25*xi.^2);
    g2 = erfc(0.5*xi) + 2*(xi.^-2).*g1;
    disp(:,i) = preFac*q*g2;
end
e_rr = (disp(3,:)-disp(1,:))./(2*dr);
e_tt = mean(disp./rvec);
e_pp = e_tt;

dStrain = [e_rr;e_tt;e_pp];
end

