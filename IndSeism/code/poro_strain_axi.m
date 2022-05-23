function dStrain = poro_strain_axi(x, y, t, q, D, Ss, b, rho_w, g, alpha, G, lambda, varargin)
% This function calculates poroelastic strain for a 2D axisymmetric case 
% using numerical solution
% INPUT:
%   r: distance from the well (m)
%   t: injection duration time (s)
%   q: injection rate (m3/s)
%   D: diffusivity (m2/s)
%   Ss: specific storage (1/m)
%   b: injection length (m)
%   rho_w: fluid density (kg/m3)
%   g: gravity acceleration (m/s2)
%   alpha: Biot's coefficient
%   G: shear modulus (Pa)
%   lambda: Lamé parameter
% OUTPUT:
%   poroelastic strain

% EQUATION: u is displacement
%               q*alpha*r*rho_w*g         /1-exp(-phi)         \ 
%       u = ----------------------------  |----------- + W(phi)|
%            8*pi*D*Ss*b*(lambda+2*G)     \    phi             /
% 
%       phi = r^2/r*D*t
% 
%       W(u) = expint(u)


% Reference
%   Rudnicki, J. W. (1986). Fluid mass sources and point forces in linear elastic diffusive solids. Mechanics of Materials, 5(4), 383–393. https://doi.org/10.1016/0167-6636(86)90042-6

%%
x = reshape(x,1,length(x));
y = reshape(y,1,length(y));

r = sqrt(x.^2+y.^2);
dr = 1e-5;
rvec = [r-dr;r;r+dr];
disp = nan(3,length(r));
K = D*alpha^2*(lambda_u+2*G)/((lambda_u-lambda)*(lambda+2*G));
preFac = q*alpha/(8*pi*K*b*(lambda+2*G));
for i = 1:length(r)
    ur = rvec(:,i).^2./(4*D*t);
    disp(:,i) = rvec(:,i)*preFac.*((1-exp(-ur))./ur+expint(ur));
end
e_rr = (disp(3,:)-disp(1,:))./(2*dr);
e_tt = mean(disp./rvec);

e_kk = e_rr+e_tt;

dStrain = [e_rr;e_tt;e_kk];
end

