function output = IndMainFunc(input,WellInfo,InjInfo,PoroElastInfo,mc,PpType,PoroType,errh,errz,signum,strike,dip,rake)
%% load parameters
lat0 = WellInfo(1);
lon0 = WellInfo(2);

q = InjInfo(:,1);
t0 = InjInfo(1,2);
dl = InjInfo(:,3);
tl = dl*86400;

D = PoroElastInfo(1);
Ss = PoroElastInfo(2);
G = PoroElastInfo(3);
nu = PoroElastInfo(4);
alpha = PoroElastInfo(5);
rho_w = PoroElastInfo(6);
g = PoroElastInfo(7);
lambda = PoroElastInfo(8);
lambda_u = PoroElastInfo(9);

K = 2*G*(1+nu)/(3*(1-2*nu));
B = alpha/(K*(Ss/(rho_w*g)));
nu_u = (3*nu+alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu));

%%
if isstruct(input)
    for istruct = 1:length(input)
        tmp = input(istruct).data;
        process;
        if istruct>1
            output = [output_old;output];
        end
        output_old = output;         
    end
else
    tmp = input;
    process;
end