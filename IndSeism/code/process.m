%% distance from the well (m)   
lat = tmp(:,3);
lon = tmp(:,2);
dep = tmp(:,4);

day = tmp(:,1)-t0;
mag = tmp(:,5);

I = day>0 & day-max(tl/86400)<30;
day = day(I);
lat = lat(I);
lon = lon(I);
mag = mag(I);
dep = dep(I);

x = distance(lat0, lon ,lat0,lon0,6378137).*sign(lon-lon0);
y = distance(lat, lon0 ,lat0,lon0,6378137).*sign(lat-lat0);

if length(errh) > 1
    errx = errh(1);
    erry = errh(2);
    errx = randn(size(lat))*errx/signum;
    erry = randn(size(lat))*erry/signum;
else
    errH = randn(size(lat))*errh/signum;
    angle2 = 2*pi*rand(size(lat));
    errx = errH.*cos(angle2);
    erry = errH.*sin(angle2);
end
errdep = randn(size(lat))*errz/signum;
x = x+errx;
y = y+erry;   
dep = dep+errdep;
d = sqrt(x.^2+y.^2);
angd = atan2(y,x);

if length(unique(mag))>1
    I = mag>=mc;
    mag = mag(I);
    d = d(I);
    angd = angd(I);
    day = day(I);
    x = x(I);
    y = y(I);
    dep = dep(I);
    
    if length(strike)>1
        strike = strike(I);
        dip = dip(I);
        rake = rake(I);
    end
end
dep0 = WellInfo(3); 
b = WellInfo(4);
dz = dep-dep0;
r = sqrt(d.^2+dz.^2);
angz = atan2(dz,d);


%% k-means earthquake density
eqdens = RadialThreeDNND_t(1,r/1e3,day,dep0/1e3,1e3);
eqdens = eqdens*10^(mc); % normalize

t = day*86400;
Nt = length(t);

t0_inj = [0;tl];
Q0_inj = [q;0];

mPp = zeros(Nt,1);
mSrr = zeros(Nt,1);
mStt = zeros(Nt,1);
mSkk = zeros(Nt,1);

for i = 1:Nt
    ItInj = t0_inj < t(i);
    vt = t(i)-t0_inj(ItInj);
    vQ = [0;Q0_inj(ItInj)];
    vQ = diff(vQ);
    Pp = 0; 
    dStrain = zeros(3,1);
    for k = 1:length(vt)
        %% Pore pressure
        if strcmp(PpType,'iso')
            Pp_tmp = Pp_iso(r(i), vt(k), vQ(k), D, Ss, rho_w, g);
        elseif strcmp(PpType,'axi')
            Pp_tmp = Pp_axi(d(i), vt(k), vQ(k), D, Ss, b, rho_w, g);
        end
        Pp = Pp+Pp_tmp';
        %% Poroelastic stress
        % Numerical, strain
        if strcmp(PoroType,'iso')
            dStrain_tmp = poro_strain_iso(r(i), vt(k), vQ(k), D, alpha, G, lambda, lambda_u);
        elseif strcmp(PoroType,'axi')
            dStrain_tmp = poro_strain_axi(d(i), 0, vt(k), vQ(k), D, b, alpha, G, lambda, lambda_u);
        end
        dStrain = dStrain+dStrain_tmp;
    end
    mPp(i) = Pp;
    
    if strcmp(PoroType,'iso')
        [Srr,Stt,Skk] = strain_to_stress_iso(dStrain,Pp, G, nu, alpha);
    elseif strcmp(PoroType,'axi')
        [Srr,Stt,Skk] = strain_to_stress_axi(dStrain,Pp, G, nu, alpha);
    end
    mSrr(i) = -Srr; % make compressional stress to be positive
    mStt(i) = -Stt;
    mSkk(i) = -Skk; 
end

taue = nan*mSrr;
for i = 1:Nt
    if length(strike)>1
        stk = strike(i);
        dp = dip(i);
        rk = rake(i);
    else
        stk = strike;
        dp = dip;
        rk = rake;
    end
    
    sigma = [mSrr(i) 0 0;0 mStt(i) 0; 0 0 mSkk(i)];
    %% local stress state to NE
    R1 = NorToGeo(pi/2-angd(i),0,-angz(i));
    SG = R1'*sigma*R1;
    
    %% NE project to Coulomb stress, coefficient of friction is 0.5
    [nn, ns, nd] = CoordVec(stk, dp, rk);
    t = SG*nn;
    td = t'*nd;
    ts = t'*ns;
    tau = sqrt(td^2+ts^2);
    sig = t'*nn;
    taue(i) = sig*0.5+tau;
end

output = table(d,r,mag,day,eqdens,mPp,taue,x,y,dep);


%% coordinates rotation
function Rng = NorToGeo(alpha, gamma, beta)
Rng = [cos(alpha)*cos(beta) sin(alpha)*cos(beta) -sin(beta);
       cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma),sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma),cos(beta)*sin(gamma);
       cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma),sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma),cos(beta)*cos(gamma)];      
end

%% fault plane vectors
function [nn, ns, nd] = CoordVec(strike, dip, rake)
nn = [-sin(strike)*sin(dip);cos(strike)*sin(dip);-cos(dip)];
ns = [cos(strike);sin(strike);0];
nd = [-sin(strike)*cos(dip);cos(strike)*cos(dip);sin(dip)];
end
