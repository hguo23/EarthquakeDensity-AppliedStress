
clear all
close all
outstr = 'output/info_Paralana.mat';

addpath('./code');
%% example: Paralana
folder = 'data';
q = 15/1e3; % injection rate
T = 176; % temperature
t0 = datenum(2011,7,10,22,9,0); % injection start time
dl = datenum(2011,7,17,12,9,0)-t0; % injection last time
lat0 = -30.213; % well latitude
lon0 = 139.72824; % well longitude
b = 200; % aquifer thickness
dep0 = 3642;

InjInfo = [q t0 dl];
WellInfo = [lat0 lon0 dep0 b];

D = 0.01; % m2/s
alpha = 0.5;
Param;
PoroElastInfo = [D, Ss, G, nu, alpha, dens, g, lambda, lambda_u];

load([folder,'/paralana.mat']);
mc = 0.2; % completeness of magnitude
PpType = 'iso'; % isotropic model
PoroType = 'iso';

errh = [33 37]; % location error (m)
errz = 48;
signum = 1; % one stardard deviation

stk = (262-180)/180*pi; % strike, from Albaric et al., 2013, Abstract, also Table 1
dip = 39/180*pi; % dip
rk = 110/180*pi; % rake

%% fitting params
conf = 90; % confidence interval

%%
Nboot = 1e1;
Nbstrp = 1e3; % bootstrapping number of fitting slope
l = 10; % bin number

for i = 1:Nboot
    output = IndMainFunc(paralana,WellInfo,InjInfo,PoroElastInfo,mc,PpType,PoroType,errh,errz,signum,stk,dip,rk);
    %%
    Eqd = output{:,5}; % density
    Pp = output{:,6}*1e6; % pore pressure
    tau = output{:,7}*1e6; % shear stress
    
    ind(i).eqd = Eqd';
    ind(i).Pp = Pp';    
    ind(i).tau = tau';
end
save(outstr,'ind');


%% display
type = 'tau'; % poroelastic stress triggered
name = 'Paralana';
load(['info_',name,'.mat']);

conf = 90;
Nbstrp = 1e3;
c = [42,156,255]/255; % color
ck = [100,100,100]/255;

    
Eqd = [ind.eqd];
Pp = [ind.Pp];
tau = [ind.tau];

eqd = Eqd;
if strcmp(type,'Pp')
    sig = Pp;
    SlopeFit_Pp;
else
    sig = tau;
    SlopeFit_Sig;
end
bin = binp;
binneg = binpneg;
binpos = binppos;
s = bvec;

[lb, hb] = ConfidenceInter(90,s);
means = median(s);

figure;hold on
errorbar(bin,p,p-p_low,p_high-p,bin-binneg,binpos-bin,'.','color',ck);
scatter(bin,p,30,'MarkerEdgeColor','k','MarkerFaceColor',c,'LineWidth',0.5);
if strcmp(type,'Pp')
    xlabel('Pore Pressure (Pa)');
else
    xlabel('Poroelastic Stress (Pa)');
end
ylabel('Earthquake Density (ev/km^3)');    
title(name);
legend(['D = ',num2str(round(D,3)),' : ',num2str(round(means,2)),' +- ',num2str(round((hb-lb)/2),3)]);
set(gca,'XScale','log','YScale','log','Fontsize',14,'TickDir','out');
box on; grid on; hold off


