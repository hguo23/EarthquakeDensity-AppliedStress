function eqdens = RadialThreeDNND_t(k,r,t,dep0,d0)
eqdens = nan(length(r),1);

for i = 1:length(r)
    It = t<t(i);
    if length(unique(r(It))) >= k        
        ktmp = k;
        dis = sort(abs(r(It)-r(i)));
        d = dis(k);
        if d == 0
            Id = dis>0;
            d = min(dis(Id));
            Id = dis == d;
            listv = 1:length(Id);
            ktmp = min(listv(Id));
        end
        if ~isempty(ktmp)
            d1 = r(i)+d;
            d2 = r(i)-d;
            V01 = 4/3*pi*d1^3;
            V02 = 4/3*pi*d2^3;            
            if d1<= dep0
                if d1<=(d0-dep0)
                    Vs1 = V01;
                else
                    H = d1-(d0-dep0);
                    Vs1 = V01-pi*H^2*(d1-H/3);
                end
            else
                if d1<(d0-dep0)
                    H = d1-dep0;
                    Vs1 = V01-pi*H^2*(d1-H/3);
                else
                    H1 = d1-dep0;
                    V1 = pi*H1^2*(d1-H1/3);
                    H2 = d1-(d0-dep0);
                    V2 = pi*H2^2*(d1-H2/3);
                    Vs1  = V01-V1-V2;
                end
            end
            if d2<= dep0
                if d2<=(d0-dep0)
                    Vs2 = V02;
                else
                    H = d2-(d0-dep0);
                    Vs2 = V02-pi*H^2*(d2-H/3);
                end
            else
                if d2<(d0-dep0)
                    H = d2-dep0;
                    Vs2 = V02-pi*H^2*(d2-H/3);
                else
                    H1 = d2-dep0;
                    V1 = pi*H1^2*(d2-H1/3);
                    H2 = d2-(d0-dep0);
                    V2 = pi*H2^2*(d2-H2/3);
                    Vs2  = V02-V1-V2;
                end
            end
            eqdens(i) = ktmp/(Vs1-Vs2);
        end
    end
end
end