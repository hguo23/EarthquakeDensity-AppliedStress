%% constants
phi = 0.2; % porosity
nu = 0.25; % poisson ratio
g = 9.81;
rho_w = 1*1e3; % fluid density
% alpha = 0.3;
G = 3e10;
K = 2*G*(1+nu)/(3*(1-2*nu));
T_K = T+273.15; 
dynVis = 1.3799566804-0.021224019151*T_K^1+1.3604562827e-4*T_K^2-4.6454090319e-7*T_K^3+8.9042735735e-10*T_K^4-9.0790692686e-13*T_K^5+3.8457331488e-16*T_K^6;
dens = 838.466135+1.40050603*T_K-0.0030112376*T_K^2+3.71822313e-7*T_K^3;

if ~exist('Cfl','var')
    Cfl = 4.2e-10;
    Cp = 1e-11;
end
Ss = (alpha/K+phi*(Cfl-Cp))*dens*g;

lambda = 2*G*nu / (1 - 2*nu);
%% temerature related

if exist('k','var') | exist('Ss','var')
    if exist('Ss','var')
        if Ss/rho_w/g < alpha/K
            error('Error. Inappropriate Input alpha or Ss.')
        end
    end    
%     if exist('k','var')
%         Ss = k*dens*g/(D*dynVis);
%         Sp = Ss/(dens*g);        
%     else
        k = dynVis*Ss*D/(dens*g);
        Sp = Ss/(dens*g);
%     end
    B = alpha/(K*Sp);
else
    error('Error. Not Enough Input Argument.')   
end
if B-1 > 1e-6
    error('Error. Inappropriate Solution of B.')   
end
Sp = alpha/(B*K);
Ss = Sp*dens*g;
k = dynVis*Ss*D/(dens*g);

nu_u = (3*nu+alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu));
lambda_u = 2*G*nu_u / (1 - 2*nu_u);
% Ku = K/(1-alpha*B);
% lambda_u1 = (3*Ku*nu_u)/(1+nu_u);
% lambda_u2 = (2*alpha*B*G+3*lambda)/3/(1 -alpha*B);
eta = alpha * (1-2*nu)/( 2*( 1 - nu));

%%
S = Ss*b; % storativity
gamma = eta/(G*S);
