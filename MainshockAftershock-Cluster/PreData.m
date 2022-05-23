clear all
close all
%%
cat = load('data/qtm_b0.8_mc1_d1.6.txt');
outfile = 'output/qtm_2D.mat';
mc = 1.0; % magnitude of completeness

t = datenum(cat(:,2),cat(:,3),cat(:,4),cat(:,5),cat(:,6),cat(:,7));
lat = cat(:,8);
lon = cat(:,9);
dep = cat(:,10);
mag = cat(:,11);
type = cat(:,12);
Id = cat(:,14);
[id, ia, ic] = unique(Id);

count = 0;
idSingle = 0;
for i = 1:length(id)
    I = ic == i;
    if sum(I)>1
        ttmp = t(I);
        lattmp = lat(I);
        lontmp = lon(I);
        deptmp = dep(I);
        magtmp = mag(I);
        typetmp = type(I);
        
        Im = typetmp == 2; % type is 2, mainshock
        Ia = typetmp == 1; % type is 1, aftershock
        Imag = magtmp>= mc;
        if sum(Ia&Imag)>0
            count = count+1;
            Rec(count).MainT = ttmp(Im);
            Rec(count).MainLat = lattmp(Im);
            Rec(count).MainLon = lontmp(Im);
            Rec(count).MainDep = deptmp(Im);
            Rec(count).MainMag = magtmp(Im);

            Rec(count).T = [ttmp(Ia&Imag)]';
            Rec(count).Lat = [lattmp(Ia&Imag)]';
            Rec(count).Lon = [lontmp(Ia&Imag)]';
            Rec(count).Dep = [deptmp(Ia&Imag)]';
            Rec(count).Mag = [magtmp(Ia&Imag)]';
        else
            idSingle = idSingle+1;
            mag_single(idSingle) = mag(id(i));
        end
    else
        idSingle = idSingle+1;
        mag_single(idSingle) = mag(id(i)); % mainshock without aftershocks
    end
end
save(outfile,'Rec','mag_single');

%% productivity analysis

rl = length(Rec);
Mag = [Rec.MainMag];
Lat = [Rec.MainLat];
Lon = [Rec.MainLon];
% Dep = [Rec.MainDep];

l = nan(1,rl);
for i = 1:rl
    l(i) = length(Rec(i).Lat); % aftershock number
end

Mag = [Mag mag_single];
l = [l mag_single*0];
magvec = min(Mag):0.1:max(Mag);
for i = 1:length(magvec)-1
    magp(i) = (magvec(i)+magvec(i+1))/2;
    I = Mag>=magvec(i) & Mag<magvec(i+1);
    lp(i) = sum(l(I))/sum(I);
end

figure;
scatter(magp,lp,50,'ro','filled')
xlabel('Magnitude');
ylabel('Number of Aftershocks');
title('Productivity Analysis');
set(gca,'Fontsize',14,'XScale','log','YScale','log');
grid on; box on;

%% find nearest focal mechanism
addpath('./code')

SelectFM;
save(outfile, 'Rec');

%% select measurement points
count = 0;
Lat = [Rec.Lat];
Lon = [Rec.Lon];
Dep = [Rec.Dep];
T = [Rec.T];

mmin = 2.5; % lower boundary of mainshock magnitude
mmax = 6; % upper boundary of mainshock magnitude

for i = 1:length(Rec)
    t = Rec(i).T;
    t = t(t>0);
    if length(t)>1 & Rec(i).MainMag<=mmax & Rec(i).MainMag>=mmin
        lat = Rec(i).Lat;
        lon = Rec(i).Lon;
        dep = Rec(i).Dep;
        mag = Rec(i).Mag;
        
        lat0 = Rec(i).MainLat;
        lon0 = Rec(i).MainLon;
        dep0 = Rec(i).MainDep;
        t0 = Rec(i).MainT;
        
        MainMag = Rec(i).MainMag;
        I = abs(Lat-lat0)<0.5 & abs(Lon-lon0)<0.5;
        latAS = Lat(I);
        lonAS = Lon(I);
        depAS = Dep(I);
        
        D = lldistkm([lat0 lon0],[latAS' lonAS']);
        D = D';
        R = sqrt(D.^2+(depAS-dep0).^2);
        
        %
        I = R<=100;
        latAS = latAS(I);
        lonAS = lonAS(I);
        depAS = depAS(I);

        if length(latAS)>1e2
            I = randi(length(latAS),1,1e2);
            latAS = latAS(I);
            lonAS = lonAS(I);
            depAS = depAS(I);
        end
        if ~isempty(Rec(i).s1)
            count = count+1;        
            Recout(count).MainMag = Rec(i).MainMag;
            Recout(count).MainLat = lat0;
            Recout(count).MainLon = lon0;
            Recout(count).MainDep = dep0;
            Recout(count).MainT = t0;
            Recout(count).Lat = latAS;
            Recout(count).Lon = lonAS;
            Recout(count).Dep = depAS;
            Recout(count).lat = lat;
            Recout(count).lon = lon;
            Recout(count).dep = dep;
            Recout(count).mag = mag;        
            Recout(count).t = t;
            Recout(count).s1 = Rec(i).s1;
            Recout(count).d1 = Rec(i).d1;
            Recout(count).r1 = Rec(i).r1;
            Recout(count).s2 = Rec(i).s2;
            Recout(count).d2 = Rec(i).d2;
            Recout(count).r2 = Rec(i).r2;
        end
    end
end
Rec = Recout;
save(outfile, 'Rec');



