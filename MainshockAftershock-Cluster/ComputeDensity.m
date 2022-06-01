clear all
close all

addpath('./codes');

%% qtm location error
errh = 0.1343*3; % in km
errz = 0.3464*3;
%%
load('output/qtm_2D.mat');
Rec = rmfield(Rec,'Lon');
Rec = rmfield(Rec,'Dep');
N  = length(Rec);

l = nan(1,N);
for i = 1:N
    l(i) = length(Rec(i).Lat);
end
Rec = rmfield(Rec,'Lat');
Rec_ori = Rec;

Nboot = 1e3; % bootstrapping number, 1000 in this study

for n = 0 :Nboot-1
    n
    %% load files
    Rec = Rec_ori;
    output = ['output/qtm_dens_',num2str(n),'.mat'];
    sta = load(['output/qtm_2D_CFS',num2str(n),'.mat']);
    
    %% load static stress output
    R1 = [sta.R1(1:l(1));sta.R1(1:l(1))];
    R2 = [sta.R2(1:l(1));sta.R2(1:l(1))];
    R1 = reshape(R1, 1, size(R1,1)*size(R1,2));
    R2 = reshape(R2, 1, size(R2,1)*size(R2,2));   
    Rec(1).R = [R1 R2];
    Rec(1).strc = [sta.stc1(1:2*l(1)) sta.stc2(1:2*l(1))];
    Rec_ori(1).errx = sta.Errx(1);
    Rec_ori(1).erry = sta.Erry(1);
    Rec_ori(1).errz = sta.Errz(1);        
    for i = 2:N
        I0 = cumsum(l(1:i-1));
        I0 = I0(end);
        Il = l(i);
        R1 = [sta.R1(I0+1:I0+Il);sta.R1(I0+1:I0+Il)];
        R2 = [sta.R2(I0+1:I0+Il);sta.R2(I0+1:I0+Il)];
        R1 = reshape(R1, 1, size(R1,1)*size(R1,2));
        R2 = reshape(R2, 1, size(R2,1)*size(R2,2));   
        Rec(i).R = [R1 R2];
        Rec(i).strc = [sta.stc1(2*I0+1:2*(I0+Il)) sta.stc2(2*I0+1:2*(I0+Il))];
        Rec_ori(i).errx = sta.Errx(i);
        Rec_ori(i).erry = sta.Erry(i);
        Rec_ori(i).errz = sta.Errz(i);        
    end
    
    %% process
    for i = 1:N        
        MainMag = Rec_ori(i).MainMag;   
        MainLat = Rec_ori(i).MainLat;
        MainLon = Rec_ori(i).MainLon;
        MainDep = Rec_ori(i).MainDep;
        lat = Rec_ori(i).lat;
        lon = Rec_ori(i).lon;
        dep = Rec_ori(i).dep;
        dt = Rec_ori(i).t - Rec_ori(i).MainT; % time to the mainshock
        
        %% distances and errors
        dx = sign(lon'-MainLon).*lldistkm([MainLat MainLon],[lat'*0+MainLat lon']);
        dy = sign(lat'-MainLat).*lldistkm([MainLat MainLon],[lat' lon'*0+MainLon]);
        errx0 = Rec_ori(i).errx;
        erry0 = Rec_ori(i).erry;
        errz0 = Rec_ori(i).errz;
        theta = 2*pi*randn(1,length(lat));
        errd = errh*randn(1,length(lat));
        errx = errd.*cos(theta);        
        erry = errd.*sin(theta);       
        d = sqrt((dx'+errx-errx0).^2+(dy'+erry-erry0).^2); 
        errdep = errz*randn(1,length(lat));        
        r = sqrt(d.^2+(dep-MainDep+errdep-errz0).^2);        
        R = Rec(i).R;

        %% stresses
        strc = Rec(i).strc;  % static stress           
        PGV = 10.^(-2.29+0.85*MainMag-1.29*log10(R)-2);        
        dyn = PGV/3.5e3*3e10; % dynamic stress
        
        %% inner and outer boundary
        r0 = 0.01*10^(0.44*MainMag)*3; % 3 times rupture length as the minimum distance
        R0 = 50; % 50 km as the maximum distance
                
        I = dt<10; % aftershocks occur in 10 days
        dt = dt(I);
        r = r(I);
        
        %% density computing
        if ~isempty(r)
            eqd = RadialThreeDNND(1,R,r,MainDep,30)/max(dt); % earthquake density
        
            I = R>r0 & R<R0 & ~isnan(eqd) & ~isinf(eqd) & eqd>0; 
            eqd = eqd(I);
            strc = strc(I);
            dyn = dyn(I);
        
            Rec(i).MainMag = MainMag;
            Rec(i).sts = strc;       
            Rec(i).dys = dyn;
            Rec(i).eqd = eqd;
        end
    end
    save(output, 'Rec');
end



