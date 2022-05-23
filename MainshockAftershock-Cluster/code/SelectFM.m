FM = load('data/sc_FM.txt');
tFM = FM(1,:);
Lat = FM(3,:);
Lon = FM(2,:);
Dep = FM(4,:);
Mag = FM(5,:);
Strike1 = FM(6,:);
Dip1 = FM(8,:);
Rake1 = FM(10,:);
Strike2 = FM(7,:);
Dip2 = FM(9,:);
Rake2 = FM(11,:);
%%
Mag = nan(1,length(Rec)*2);
strike = nan(1,length(Rec)*2);
dip = nan(1,length(Rec)*2);
rake = nan(1,length(Rec)*2);

for i = 1:length(Rec)
    lat0 = Rec(i).MainLat;
    lon0 = Rec(i).MainLon;
    mag = Rec(i).MainMag;
    dis = lldistkm([lat0 lon0],[Lat' Lon']); % epicentral distance
    dep0 = Rec(i).MainDep;
    dis = sqrt(dis.^2+(dep0-Dep').^2); % hypocentral distance
    [v,Il] = min(dis); % find the nearest focal mechanism solution
   
    if v<10 % if the solution is less than 10 km
        strike(2*i-1:2*i) = [Strike1(Il) Strike2(Il)];
        dip(2*i-1:2*i) = [Dip1(Il) Dip2(Il)];
        rake(2*i-1:2*i) = [Rake1(Il) Rake2(Il)];
        Mag(2*i-1:2*i) = [mag mag];
    end
end

I = dip == 90;
dip(I) = 89;
I = sind(strike) == 0;
strike(I) = strike(I)-0.1;

strike = reshape(strike,2,length(Rec));
dip = reshape(dip,2,length(Rec));
rake = reshape(rake,2,length(Rec));

for i = 1:length(Rec)
    if ~isnan(strike(1,i)) & ~isnan(strike(2,i))
        Rec(i).s1 = strike(1,i);
        Rec(i).d1 = dip(1,i);
        Rec(i).r1 = rake(1,i);
        Rec(i).s2 = strike(2,i);
        Rec(i).d2 = dip(2,i);
        Rec(i).r2 = rake(2,i);
    end
end