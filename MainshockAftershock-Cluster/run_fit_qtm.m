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
    eqd = [Rec.eqd]; % radial density
   
    I = sts>0 ;
    sts = sts(I);
    eqd2 = eqd(I);
        
    slope(n+1).sts = sts;
    slope(n+1).dys = dys;
    slope(n+1).eqd1 = eqd;
    slope(n+1).eqd2 = eqd2;    
end

%%
clc
clearvars -except slope
sts = [slope.sts];
dys = [slope.dys]*7.07945784384137^4;
eqd1 = [slope.eqd1];
eqd2 = [slope.eqd2];

n = 16;
ld = max([round(length(dys)/n),1e3]); 
sld = ld;
ls = max([round(length(sts)/n),1e3]); 
sls = ls;
[bindp,bindneg,bindpos,p1, p1_low, p1_high,dfe1] = calP_window(dys,eqd1,ld,sld);
[binsp,binsneg,binspos,p2, p2_low, p2_high,dfe2] = calP_window(sts,eqd2,ls,sls);
     
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

%% slope values in 90% confidence interval
[LB1, UB1] = ConfidenceInter(90,bvec1);
[LB2, UB2] = ConfidenceInter(90,bvec2);
