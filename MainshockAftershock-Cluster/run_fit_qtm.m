clear all
close all
addpath('./code');
%%
Nboot = 1; % bootstrapping number in analyzing location errors
name = 'QTM';
%% 
for n = 0:Nboot-1
    i = n
    load(['output/qtm_dens_',num2str(i),'.mat']);
        
    dys = [Rec.dys]; % normalized dynamic stresses
    sts = [Rec.sts]; % static stresses    
    eqd = [Rec.eqd]; % radial density for static stresses analysis
    eqdn = [Rec.eqdn]; % radial density for dynamic stresses analysis
    
    I = sts>0 & eqd>1e-6;
    sts = sts(I);
    eqd = eqd(I);
    
    I = dys>0 & eqdn>1e-8;
    dys = dys(I);
    eqdn = eqdn(I);
    
    slope(n+1).sts = sts;
    slope(n+1).dys = dys;
    slope(n+1).eqd = eqd;
    slope(n+1).eqdn = eqdn;    
end

%%
clc
clearvars -except slope
sts = [slope.sts];
dys = [slope.dys]*7.07945784384137^4;
eqd = [slope.eqd];
eqdn = [slope.eqdn];

n = 25;
ld = max([round(length(dys)/n),1e3]); 
sld = ld;
ls = max([round(length(sts)/n),1e3]); 
sls = ls;
[bindp,bindneg,bindpos,p1, p1_low, p1_high,dfe1] = calP_window(dys,eqdn,ld,sld);
[binsp,binsneg,binspos,p2, p2_low, p2_high,dfe2] = calP_window(sts,eqd,ls,sls);
     
%% fit slope
conf = 90; % confidence interval
Nbstrp = 1e3; % bootstrapping number in fitting slope
c1 = [[3,67,128]/255;[127,176,240]/255]; % blue color
c2 = [[128,11,3]/255;[235,158,52]/255]; % red color
 
if length(p1(p1>0)) > 3 % at least 3 bins are not empty
    [~,~,~,~,bvec1] = FitSlope_error(bindp,p1,p1_low,p1_high,dfe1,conf,Nbstrp);
else
    bvec1 = [];
end
if length(p2(p2>0)) > 3 % at least 3 bins are not empty
    [~,~,~,~,bvec2] = FitSlope_error(binsp,p2,p2_low,p2_high,dfe2,conf,Nbstrp);
else
    bvec2 = [];
end

% display binning result    
figure(2);clf; hold on
errorbar(binsp,p2,p2-p2_low,p2_high-p2,binsp-binsneg,binspos-binsp,'.','color',c1(1,:),'LineWidth',1,'MarkerEdgeColor',c1(1,:),'MarkerFaceColor',c1(1,:));
errorbar(bindp,p1,p1-p1_low,p1_high-p1,bindp-bindneg,bindpos-bindp,'.','color',c2(1,:),'LineWidth',1,'MarkerEdgeColor',c2(1,:),'MarkerFaceColor',c2(1,:));
h1 = scatter(bindp,p1,30,c2(1,:),'s','filled');
h2 = scatter(binsp,p2,30,c1(1,:),'s','filled');
xlabel('Stress (Pa)');
ylabel('Normalized Earthquake Density (ev/km^3)');
legend([h1 h2],{'dynamic stress','static stress'},'location','southeast');
set(gca,'XScale','log','YScale','log','Fontsize',14);
box on; grid on; hold off

figure(3);hist(bvec1)
figure(4);hist(bvec2)
